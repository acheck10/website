if (!exists(".src_cpp")) source("R/rcpp_cache.R")      # persistent, self-repairing compile cache for the C++ kernels
if (!exists("ham_ll")) .src_cpp("R/hamilton.cpp")      # C++ kernels: Hamilton filter, ffbs_cpp, ffbs_rw_cpp, sv_ind_cpp, dens_norm2, dens_t2
# Gibbs samplers for two-state Markov-switching models (Piger 2026)
#   msm_gibbs(): switching mean,     y_t = a_{S_t} + sigma * eps_t,  eps ~ N(0,1)
#   msq_gibbs(): switching quantile, y_t = a_{S_t} + sigma * w_t,    w ~ AL(0,1,tau)
# State 1 = low regime (a1 < a0). Both use the same priors (Piger Table 1).
# Base R only.

# ---- helpers ---------------------------------------------------------------

# Inverse-Gaussian draws (Michael, Schucany & Haas 1976); mu, lambda vectorised
rinvgauss_vec <- function(n, mu, lambda) {
  nu <- rnorm(n)^2
  x  <- mu + mu^2 * nu / (2 * lambda) -
        mu / (2 * lambda) * sqrt(4 * mu * lambda * nu + mu^2 * nu^2)
  u  <- runif(n)
  ifelse(u <= mu / (mu + x), x, mu^2 / x)
}

# Forward-filter / backward-sample for a 2-state chain.
# dens: T x 2 matrix of f(y_t | S_t = j), column 1 -> S=0, column 2 -> S=1.
# P00 = Pr(S_t=0|S_{t-1}=0), P11 = Pr(S_t=1|S_{t-1}=1).
# Returns list(S = draw of states in {0,1}, filt = filtered Pr(S_t=1|Y_t)).
ffbs_R <- function(dens, P00, P11, init1 = NULL) {      # reference implementation (the C++ ffbs_cpp below reproduces it)
  n <- nrow(dens)
  Tm <- matrix(c(P00, 1 - P00, 1 - P11, P11), 2, 2, byrow = TRUE)  # row=from
  # initial distribution: ergodic by default, or a fixed Pr(S_1 = 1) = init1
  pi1 <- if (is.null(init1)) (1 - P00) / ((1 - P00) + (1 - P11)) else init1
  pred <- c(1 - pi1, pi1)
  filt <- matrix(0, n, 2)
  for (t in 1:n) {
    if (t > 1) pred <- as.numeric(filt[t - 1, ] %*% Tm)
    jt <- pred * dens[t, ]
    s  <- sum(jt)
    if (!is.finite(s) || s <= 0) jt <- pred else jt <- jt / s
    filt[t, ] <- jt
  }
  S <- integer(n)
  S[n] <- as.integer(runif(1) < filt[n, 2])
  for (t in (n - 1):1) {
    nxt <- S[t + 1] + 1
    w0  <- filt[t, 1] * Tm[1, nxt]
    w1  <- filt[t, 2] * Tm[2, nxt]
    S[t] <- as.integer(runif(1) < w1 / (w0 + w1))
  }
  list(S = S, filt = filt[, 2])
}

ffbs <- function(dens, P00, P11, init1 = NULL) ffbs_cpp(dens, P00, P11, if (is.null(init1)) -1 else init1)

# Draw (a0, a1) from N(Vbar * (sum W_t y_t w_t + V^-1 mu), Vbar), w_t = precision weight
draw_alpha <- function(y, S, w, mu, Vinv, a1max = NULL) {
  # a1max: optional upper bound on the low-regime mean a1 (e.g. 0 = 'recession regime must have negative growth').
  # W_t = (1-S_t, S_t)'; with two disjoint regimes the matrix is diagonal
  p0 <- sum(w[S == 0]);  p1 <- sum(w[S == 1])
  q0 <- sum((w * y)[S == 0]); q1 <- sum((w * y)[S == 1])
  prec <- matrix(c(p0, 0, 0, p1), 2, 2) + Vinv
  Vbar <- solve(prec)
  mbar <- Vbar %*% (c(q0, q1) + Vinv %*% mu)
  L <- t(chol(Vbar))
  # Regime normalisation a1 < a0 (and a1 < a1max if given) as a truncation of the normal conditional.
  ub <- if (is.null(a1max)) Inf else a1max
  if (abs(Vbar[1, 2]) < 1e-12) {
    # Independent normals (diagonal prior and data precision): exact two-step draw that works even far in the tail.
    #   a1 ~ N(m1, s1^2) truncated to (-Inf, ub), accepted with probability P(a0 > a1); then a0 | a1 ~ N(m0, s0^2) truncated to (a1, Inf).
    m0 <- mbar[1]; s0 <- sqrt(Vbar[1, 1]); m1 <- mbar[2]; s1 <- sqrt(Vbar[2, 2])
    lcdf_ub <- pnorm((ub - m1) / s1, log.p = TRUE)
    for (k in 1:100000) {
      a1 <- if (is.infinite(ub)) m1 + s1 * rnorm(1) else m1 + s1 * qnorm(log(runif(1)) + lcdf_ub, log.p = TRUE)
      lsurv <- pnorm((a1 - m0) / s0, lower.tail = FALSE, log.p = TRUE)       # log P(a0 > a1)
      if (log(runif(1)) < lsurv) {
        a0 <- m0 + s0 * qnorm(log(runif(1)) + lsurv, lower.tail = FALSE, log.p = TRUE)
        return(c(a0, a1))
      }
    }
  }
  for (k in 1:100000) { a <- as.numeric(mbar + L %*% rnorm(2)); if (a[2] < a[1] && a[2] < ub) return(a) }
  a[2] <- min(a[1], ub) - abs(rnorm(1, 0, 0.01)); a      # last-resort fallback (non-diagonal prior only)
}

draw_P <- function(S, prior) {
  n <- length(S)
  N00 <- sum(S[-n] == 0 & S[-1] == 0); N01 <- sum(S[-n] == 0 & S[-1] == 1)
  N10 <- sum(S[-n] == 1 & S[-1] == 0); N11 <- sum(S[-n] == 1 & S[-1] == 1)
  c(P00 = rbeta(1, N00 + prior$B00, N01 + prior$B01),
    P11 = rbeta(1, max(N11 + prior$B11, 1e-6), N10 + prior$B10))
}

# Transition-probability priors by data frequency. QUARTERLY (Eo & Kim): P00 ~ Beta(9, 1), P11 ~ Beta(8, 2): mean regime durations 10 and 5 quarters, 10 pseudo-observations each.
# MONTHLY: the same durations in CALENDAR time (30 and 15 months) at the same strength relative to the sample size (three times as many months per year, so three times the pseudo-observations): P00 ~ Beta(29, 1), P11 ~ Beta(28, 2).
trans_prior <- function(freq = c("quarterly", "monthly")) if (match.arg(freq) == "monthly") list(B00 = 29, B01 = 1, B11 = 28, B10 = 2) else list(B00 = 9, B01 = 1, B11 = 8, B10 = 2)

# Regime means: a0~N(.5,1), a1~N(-.5,1) (Piger Table 1). Transition probabilities: P00~Beta(9,1), P11~Beta(8,2),
# the same proper priors Eo & Kim use for q and p (means .9 and .8, about 10 pseudo-observations each).
# NOTE: Table 1 of Piger lists P11 ~ Beta(0.0, 1.0). Read as (stay, leave) pseudo-counts that is Beta(0, 1), which
# is degenerate at 0: on the GDP data it let the low regime become absorbing at P11 = 0 (one-quarter "recessions")
# for MSM-t, the SV models and MSQ(0.5). The earlier default (B00 = 9.5, B01 = .5, B11 = 0, B10 = 1) is kept
# in the comment so results can be traced: use modifyList(default_prior, list(B00 = 9.5, B01 = .5, B11 = 0, B10 = 1)).
default_prior <- list(mu = c(0.5, -0.5), V = diag(c(1, 1)),
                      B00 = 9, B01 = 1, B11 = 8, B10 = 2,
                      sig_a = 0, sig_b = 0)   # sigma prior IG(a,b); a=b=0 is the improper 1/sigma^2 (or 1/sigma) default
# Samplers also accept `start` (a state list to resume from; used by the Geweke test)
# and `guard` (TRUE = avoid empty regimes; set FALSE for exact-posterior tests).
# Each returns its final state in $state.

# ---- Markov-switching mean -------------------------------------------------

msm_gibbs <- function(y, ndraw = 5000, burn = 1000, prior = default_prior,
                      seed = NULL, verbose = TRUE, start = NULL, guard = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y); Vinv <- solve(prior$V)
  if (is.null(start)) {
    S <- as.integer(y < quantile(y, 0.15))
    a <- c(mean(y[S == 0]), mean(y[S == 1])); sig2 <- var(y)
    P <- c(P00 = 0.95, P11 = 0.8)
  } else { S <- start$S; a <- start$a; sig2 <- start$sig2; P <- start$P }
  keep <- ndraw - burn
  out <- list(alpha = matrix(NA, keep, 2), sigma = numeric(keep),
              P = matrix(NA, keep, 2), Sprob = numeric(n))
  for (i in 1:ndraw) {
    # states
    dens <- cbind(dnorm(y, a[1], sqrt(sig2)), dnorm(y, a[2], sqrt(sig2)))
    S <- ffbs(dens, P[1], P[2], prior$init_p1)$S
    if (guard && (all(S == 0) || all(S == 1))) S[which.min(y)] <- 1L - S[1]  # avoid empty regime
    # alpha | S, sigma^2
    a <- draw_alpha(y, S, rep(1 / sig2, n), prior$mu, Vinv, prior$a1_max)
    # sigma^2 | alpha, S : IG(T/2, SSR/2)  (prior 1/sigma^2)
    ssr  <- sum((y - a[S + 1])^2)
    sig2 <- 1 / rgamma(1, n / 2 + prior$sig_a, rate = ssr / 2 + prior$sig_b)
    P <- draw_P(S, prior)
    # label normalisation a1 < a0
    if (a[2] > a[1]) { a <- rev(a); S <- 1L - S; P <- c(P[2], P[1]) }
    if (i > burn) {
      k <- i - burn
      out$alpha[k, ] <- a; out$sigma[k] <- sqrt(sig2); out$P[k, ] <- P
      out$Sprob <- out$Sprob + S / keep
    }
    if (verbose && i %% 1000 == 0) message("MSM draw ", i, "/", ndraw)
  }
  out$state <- list(S = S, a = a, sig2 = sig2, P = P)
  out
}

# ---- Markov-switching quantile ---------------------------------------------

msq_gibbs <- function(y, tau = 0.3, ndraw = 5000, burn = 1000,
                      prior = default_prior, seed = NULL, verbose = TRUE,
                      start = NULL, guard = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y); Vinv <- solve(prior$V)
  theta <- (1 - 2 * tau) / (tau * (1 - tau))
  psi2  <- 2 / (tau * (1 - tau))
  if (is.null(start)) {
    S <- as.integer(y < quantile(y, 0.15))
    a <- c(quantile(y[S == 0], tau), quantile(y[S == 1], tau)); names(a) <- NULL
    sig <- sd(y); zs <- rep(sig, n)         # z* ~ Exp(mean sigma)
    P <- c(P00 = 0.95, P11 = 0.8)
  } else { S <- start$S; a <- start$a; sig <- start$sig; zs <- start$zs; P <- start$P }
  keep <- ndraw - burn
  out <- list(alpha = matrix(NA, keep, 2), sigma = numeric(keep),
              P = matrix(NA, keep, 2), Sprob = numeric(n))
  for (i in 1:ndraw) {
    # states | z*: y_t|S,z* ~ N(a_S + theta z*, psi2 * sigma * z*)
    sd_t <- sqrt(psi2 * sig * zs)
    dens <- cbind(dnorm(y, a[1] + theta * zs, sd_t), dnorm(y, a[2] + theta * zs, sd_t))
    S <- ffbs(dens, P[1], P[2], prior$init_p1)$S
    if (guard && (all(S == 0) || all(S == 1))) S[which.min(y)] <- 1L - S[1]
    # alpha | S, z*, sigma
    a <- draw_alpha(y - theta * zs, S, 1 / (psi2 * sig * zs), prior$mu, Vinv, prior$a1_max) |> as.numeric()
    # sigma | ... : IG(3T/2, sum[(y-a-theta z*)^2/(2 psi2 z*) + z*])  (prior 1/sigma)
    r   <- y - a[S + 1]
    rate <- sum((r - theta * zs)^2 / (2 * psi2 * zs) + zs)
    sig <- 1 / rgamma(1, 3 * n / 2 + prior$sig_a, rate = rate + prior$sig_b)
    # z* | ... : GIG(1/2, A, B); 1/z* ~ InvGauss(sqrt(A/B), A)
    A <- theta^2 / (psi2 * sig) + 2 / sig
    B <- pmax(r^2 / (psi2 * sig), 1e-10)
    zs <- 1 / rinvgauss_vec(n, sqrt(A / B), A)
    P <- draw_P(S, prior)
    if (a[2] > a[1]) { a <- rev(a); S <- 1L - S; P <- c(P[2], P[1]) }
    if (i > burn) {
      k <- i - burn
      out$alpha[k, ] <- a; out$sigma[k] <- sig; out$P[k, ] <- P
      out$Sprob <- out$Sprob + S / keep
    }
    if (verbose && i %% 1000 == 0) message("MSQ(", tau, ") draw ", i, "/", ndraw)
  }
  out$state <- list(S = S, a = a, sig = sig, zs = zs, P = P)
  out
}

# ---- classification metrics (Piger Sec. 4) -----------------------------------
# prob: classifier in [0,1]; nber: 0/1 recession indicator. Threshold = Youden max.
classify_metrics <- function(prob, nber) {
  cs  <- sort(unique(c(0, prob, 1)))
  tpr <- sapply(cs, function(c) mean(prob[nber == 1] >= c))
  fpr <- sapply(cs, function(c) mean(prob[nber == 0] >= c))
  o <- order(fpr, tpr)
  auc <- sum(diff(fpr[o]) * (head(tpr[o], -1) + tail(tpr[o], -1)) / 2)
  cstar <- cs[which.max(tpr - fpr)]
  pred <- prob >= cstar
  prec <- if (sum(pred) > 0) mean(nber[pred] == 1) else NA
  tp <- mean(prob[nber == 1] >= cstar)
  c(AUC = auc, TPR = tp, FPR = mean(prob[nber == 0] >= cstar), Precision = prec,
    F1 = 2 * tp * prec / (tp + prec))
}

# Prior on the stochastic-volatility innovation variance sigma2_eta ~ IG(nu_eta/2, h_eta/2). The paper's printed IG(4,1) (nu = 8, h = 2) is incompatible with its reported
# posterior (and with the prior its own code uses, IG(5,.05)); the DEFAULT here is a nearly flat IG(1,.01) (nu = 2, h = .02). Set the environment variable SV_PRIOR=paper for IG(4,1).
sv_eta_prior_default <- if (Sys.getenv("SV_PRIOR", "flat") == "paper") c(8, 2) else c(2, 0.02)
