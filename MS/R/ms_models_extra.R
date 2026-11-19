if (!exists(".src_cpp")) source("R/rcpp_cache.R")      # persistent, self-repairing compile cache for the C++ kernels
if (!exists("ham_ll")) .src_cpp("R/hamilton.cpp")
# Two additional candidate models. Requires ms_samplers.R (ffbs, draw_alpha, draw_P, ...).
#
#  msmt_gibbs()   : Markov-switching mean with Student-t errors
#                     y_t = a_{S_t} + sigma * e_t,  e_t ~ t_nu
#                   via the normal scale mixture e_t = z_t / sqrt(lambda_t), lambda_t ~ Gamma(nu/2, nu/2).
#                   nu is either fixed or drawn on a grid (uniform prior).
#
#  eokim_gibbs()  : Eo & Kim (2013), "Markov-Switching Models with Evolving Regime-Specific
#                   Parameters". Their Table 3 model:
#                     y_t     = delta_t + (1-S_t) mu0*_t + S_t mu1*_t + e_t,  e_t ~ N(0, exp(h_t))
#                     mu_j*_t = mu_j*_{t-1} + theta_j * 1{switch at t} * (pi0 mu0*_{t-1} + pi1 mu1*_{t-1}) + w_jt
#                               w_0t ~ N(0, d10_t s2w0),  w_1t ~ N(0, d01_t s2w1)
#                     delta_t = delta_{t-1} + eps_t,  Var(eps) fixed at 0.0004
#                     h_t     = h_{t-1} + eta_t,      eta ~ N(0, s2eta)
#                     mu0* > 0, mu1* < 0, -1 < 1 + theta0 pi0 + theta1 pi1 < 1
#                   Priors as in Eo-Kim Table 3 (moments matched to inverse-gamma / beta).
#
# Deviations from the paper's algorithm (see notes in the page / README):
#   * Step 1 samples mu* with a Carter-Kohn smoother whose singular state noise is regularised by a
#     tiny ridge (kappa), instead of the paper's special-case rules for switch dates.
#   * Stochastic volatility is sampled with the Kim-Shephard-Chib 7-component mixture.
#   * The paper's initial priors for mu*_0 and delta_0 are not reported; N(.,1) priors are used.

# ---------------------------------------------------------------------------
# MSM with t errors
# ---------------------------------------------------------------------------
msmt_gibbs <- function(y, ndraw = 5000, burn = 1000, prior = default_prior,
                       est_nu = FALSE, nu = 5, nu_grid = c(3, 4, 5, 6, 8, 10, 15, 20, 30),
                       seed = NULL, verbose = TRUE, start = NULL, guard = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y); Vinv <- solve(prior$V)
  if (is.null(start)) {
    S <- as.integer(y < quantile(y, 0.15))
    a <- c(mean(y[S == 0]), mean(y[S == 1])); sig2 <- var(y)
    lam <- rep(1, n); P <- c(P00 = 0.95, P11 = 0.8)
  } else { S <- start$S; a <- start$a; sig2 <- start$sig2; lam <- start$lam; if (est_nu) nu <- start$nu; P <- start$P }
  keep <- ndraw - burn
  out <- list(alpha = matrix(NA, keep, 2), sigma = numeric(keep), nu = numeric(keep),
              P = matrix(NA, keep, 2), Sprob = numeric(n))
  for (i in 1:ndraw) {
    sig <- sqrt(sig2)
    dens <- cbind(dt((y - a[1]) / sig, nu) / sig, dt((y - a[2]) / sig, nu) / sig)
    S <- ffbs(dens, P[1], P[2], prior$init_p1)$S
    if (guard && (all(S == 0) || all(S == 1))) S[which.min(y)] <- 1L - S[1]
    # S was drawn with lambda integrated out, so lambda must be redrawn given the new S
    # (partially collapsed Gibbs: (S, lambda) jointly, then the parameters)
    lam  <- rgamma(n, (nu + 1) / 2, rate = (nu + (y - a[S + 1])^2 / sig2) / 2)
    a <- draw_alpha(y, S, lam / sig2, prior$mu, Vinv, prior$a1_max)
    r <- y - a[S + 1]
    sig2 <- 1 / rgamma(1, n / 2 + prior$sig_a, rate = sum(lam * r^2) / 2 + prior$sig_b)
    if (est_nu) {
      lp <- sapply(nu_grid, function(v) sum(dgamma(lam, v / 2, rate = v / 2, log = TRUE)))
      nu <- sample(nu_grid, 1, prob = exp(lp - max(lp)))
    }
    P <- draw_P(S, prior)
    if (a[2] > a[1]) { a <- rev(a); S <- 1L - S; P <- c(P[2], P[1]) }
    if (i > burn) {
      k <- i - burn
      out$alpha[k, ] <- a; out$sigma[k] <- sqrt(sig2); out$nu[k] <- nu; out$P[k, ] <- P
      out$Sprob <- out$Sprob + S / keep
    }
    if (verbose && i %% 1000 == 0) message("MSM-t draw ", i, "/", ndraw)
  }
  out$state <- list(S = S, a = a, sig2 = sig2, lam = lam, nu = nu, P = P)
  out
}

# ---------------------------------------------------------------------------
# Eo & Kim (2013)
# ---------------------------------------------------------------------------
eokim_prior <- list(theta_m = -0.1, theta_v = 0.25,      # theta_j ~ N(-0.1, 0.5^2)
                    nu_w = 10, h_w = 2,                   # s2w ~ IG(5,1): mean .25, sd .144
                    nu_eta = sv_eta_prior_default[1], h_eta = sv_eta_prior_default[2],   # s2eta: default nearly flat IG(1,.01) (SV_PRIOR=paper: the printed IG(4,1))
                    B00 = 9, B01 = 1, B11 = 8, B10 = 2)   # q~Beta(9,1), p~Beta(8,2)

# KSC (1998) 7-component normal mixture approximation to log chi^2_1
.ksc_m <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) - 1.2704
.ksc_v <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)
.ksc_q <- c(0.0073, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)

# Omori, Chib, Shephard & Nakajima (2007) 10-component mixture for log chi^2_1
# (columns: weight, mean, variance; means already include the -1.2704 shift).
# Same table as in carter_kohn_SV_draw.jl.
.omori10 <- cbind(
  p = c(0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575, 0.00115),
  m = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65),
  v = c(0.11265, 0.17788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342))
.ksc7 <- cbind(p = .ksc_q, m = .ksc_m, v = .ksc_v)

# One stochastic-volatility draw: e_t = exp(h_t/2) * N(0,1) is the residual, h_t a random walk.
# Draws mixture indicators for y*_t = log(e_t^2 + cc) - h_t, then h | indicators (Carter-Kohn).
# bug: test-only variants that confuse a standard deviation with a variance in the log-chi^2 mixture step (mix[, 3] is a VARIANCE):
#   "ind_v"  : indicator probabilities use mix[, 3] as the standard deviation;  "kal_sqrt": the Kalman smoother uses sqrt(mix[, 3]) as the noise variance;
#   "kal_sq" : the Kalman smoother uses mix[, 3]^2 as the noise variance;       "ind_v_kal_sq": both ind_v and kal_sq (code treats the column as an SD).
sv_draw <- function(e, h, s2eta, h0m, h0v, mix = .omori10, cc = 1e-8, bug = "") {
  ys <- log(e^2 + cc); n <- length(ys); K <- nrow(mix)
  if (bug == "") { kk <- sv_ind_cpp(ys, h, mix); return(ffbs_rw(ys - mix[kk, 2], mix[kk, 3], s2eta, h0m, h0v)) }      # fast path (C++ indicator draw); the R code below is kept for the test-only 'bug' variants
  pm <- sapply(1:K, function(k) mix[k, 1] * dnorm(ys - h, mix[k, 2], if (bug %in% c("ind_v", "ind_v_kal_sq")) mix[k, 3] else sqrt(mix[k, 3])))
  pm <- matrix(pm, n, K); rs <- rowSums(pm); bad <- !is.finite(rs) | rs <= 0
  if (any(bad)) pm[bad, ] <- rep(mix[, 1], each = sum(bad))
  pm <- pm / rowSums(pm)
  kk <- pmin(1 + rowSums(runif(n) > t(apply(pm, 1, cumsum))), K)
  Rk <- if (grepl("^kal_x", bug)) mix[kk, 3] * as.numeric(sub("kal_x", "", bug)) else if (bug == "kal_sqrt") sqrt(mix[kk, 3]) else if (bug %in% c("kal_sq", "ind_v_kal_sq")) mix[kk, 3]^2 else mix[kk, 3]
  ffbs_rw(ys - mix[kk, 2], Rk, s2eta, h0m, h0v)
}


# Single-move Metropolis update of the log-volatility path (Kim, Shephard & Chib 1998, Sec. 2.2-style): h_t is updated one period at a time given its neighbours,
# with the proposal built from a first-order expansion of the likelihood. This mixes much more slowly than the block (mixture + Carter-Kohn) draw in sv_draw();
# used only to test whether a slow single-move sampler could explain the paper's smooth volatility path and small s2eta.
sv_single_move <- function(e, h, s2eta, h0m, h0v) {
  n <- length(h); e2 <- e^2
  for (t in 1:n) {
    if (t == 1) { pv <- 1 / (1 / h0v + 1 / s2eta); mu <- pv * (h0m / h0v + h[2] / s2eta) }
    else if (t == n) { pv <- s2eta; mu <- h[n - 1] }
    else { pv <- s2eta / 2; mu <- (h[t - 1] + h[t + 1]) / 2 }
    mt <- mu + (pv / 2) * (e2[t] * exp(-mu) - 1)                         # proposal mean (linearised likelihood)
    hn <- rnorm(1, mt, sqrt(pv))
    lw <- function(x) (-x / 2 - e2[t] * exp(-x) / 2) - (-x / 2 - e2[t] * exp(-mu) * (1 - (x - mu)) / 2)   # log(true / linearised likelihood)
    if (log(runif(1)) < lw(hn) - lw(h[t])) h[t] <- hn
  }
  h
}

# scalar random-walk-plus-noise FFBS: z_t = h_t + N(0, R_t), h_t = h_{t-1} + N(0, s2eta)
ffbs_rw_R <- function(z, R, s2eta, h0m, h0v) {      # reference implementation (ffbs_rw_cpp reproduces it)
  n <- length(z); m <- numeric(n); P <- numeric(n)
  mp <- h0m; Pp <- h0v
  for (t in 1:n) {
    if (t > 1) { mp <- m[t - 1]; Pp <- P[t - 1] + s2eta }
    K <- Pp / (Pp + R[t]); m[t] <- mp + K * (z[t] - mp); P[t] <- (1 - K) * Pp
  }
  h <- numeric(n); h[n] <- rnorm(1, m[n], sqrt(P[n]))
  for (t in (n - 1):1) {
    G <- P[t] / (P[t] + s2eta)
    h[t] <- rnorm(1, m[t] + G * (h[t + 1] - m[t]), sqrt(max(P[t] * (1 - G), 1e-12)))
  }
  h
}
ffbs_rw <- function(z, R, s2eta, h0m, h0v) ffbs_rw_cpp(z, if (length(R) == 1) rep(R, length(z)) else R, s2eta, h0m, h0v)

eokim_gibbs <- function(y, ndraw = 3000, burn = 1000, prior = eokim_prior,
                        sig2_eps = 0.0004, init_mu = c(0.5, -0.5), kappa = 1e-9,
                        seed = NULL, verbose = TRUE, start = NULL, guard = TRUE,
                        delta0 = NULL, h0m = NULL, h0v = 5, init_p1 = NULL,
                        signs = TRUE, qp_constraint = TRUE,
                        update_S = TRUE, update_qp = TRUE, update_theta = TRUE, update_h = TRUE,
                        sv_mix = .omori10, sv_c = 1e-8, init_q = 0.15, delta_sd = 1, mu_sd = 1, paper_rules = FALSE, const_var = FALSE, sv_sd_mode = FALSE, sv_prior_sd = FALSE, sv_single = FALSE, sv_bug = "", s2eta_incr = "log", use_cpp = TRUE, t_errors = FALSE, nu = 5) {
  # t_errors = TRUE: Student-t errors with nu degrees of freedom (fixed) through the scale mixture e_t / sqrt(lam_t), lam_t ~ Gamma(nu/2, nu/2): the measurement variance becomes exp(h_t) / lam_t in every step, lam_t | S, X, h is a Gamma draw.
  # signs = FALSE additionally relabels each retained draw into the identified labelling (mean(mu0*) > mean(mu1*)): without the sign restrictions the regime labels are not identified.
  # sig2_eps = 0: constant long-run mean (no random-walk drift).
  # use_cpp: Carter-Kohn draw and single-move regime update in C++ (eokim_ck_cpp / eokim_S_cpp in hamilton.cpp, same algorithm and random-number order; FALSE = the R loops below). s2eta_incr: which increments enter the sigma2_eta update: "log" (correct: diff of ln sigma2_t), "level" (diff of sigma2_t), "sd" (diff of sigma_t). Test-only.
  # sv_single = TRUE: update h by the single-move Metropolis sampler sv_single_move() instead of the block mixture sampler.
  # sv_prior_sd = TRUE: the IG(nu_eta/2, h_eta/2) prior is placed on the innovation STANDARD DEVIATION sigma_eta (not its variance); sigma_eta is drawn from its
  # (non-conjugate) full conditional on a grid and s2eta = sigma_eta^2 is what enters the state equation.
  # sv_sd_mode = TRUE: use the drawn s2eta as the random-walk innovation STANDARD DEVIATION in the volatility state equation (variance s2eta^2) while still
  # updating it with the usual inverse-gamma step. A test of whether the paper's sampler mixes up sd and variance (investigation/FINDINGS_EK_gap.md).
  # mu_sd: prior sd of the initial mu0*, mu1*. paper_rules: use Eo-Kim's literal backward-sampling rules (i)-(iii) (not an exact
  # draw: the active/other component at switch dates is set to its conditional mean). const_var: constant error variance instead of SV.
  # delta_sd: prior sd of the initial long-run growth delta_1 ~ N(delta0, delta_sd^2) (default 1 = diffuse).
  # init_q: the initial S marks the lowest init_q share of y as recession (used to start chains from different points).
  # sv_mix / sv_c: mixture table and offset for the stochastic-volatility draw. The default
  # (Omori et al. 10-component, offset 1e-8) passes a Geweke test; the earlier KSC 7-component
  # table with offset 1e-3 does not (see geweke_sv.R).
  # update_S / update_qp / update_theta = FALSE hold that block fixed (used only to localise
  # problems with the Geweke test; the production sampler updates everything).
  # signs: impose mu0* > 0, mu1* < 0 (paper eq. 35). qp_constraint: also impose the stationarity
  # restriction (38) when drawing (q, p), since it involves pi0, pi1 (theta draws always impose it).
  # delta0 / h0m: prior means for delta_1 and h_1. NULL = data-based (ybar, log var(y));
  # the Geweke test passes constants because priors must not depend on y.
  if (!is.null(seed)) set.seed(seed)
  n <- length(y); ybar <- mean(y)
  if (is.null(delta0)) delta0 <- ybar
  if (is.null(h0m)) h0m <- log(var(y))
  # ---- initial values --------------------------------------------------
  if (is.null(start)) {
    S <- as.integer(y < quantile(y, init_q))
    theta <- c(-0.1, -0.1); s2w <- c(0.25, 0.25); s2eta <- 0.05
    q <- 0.9; p <- 0.8
    h <- rep(log(var(y)), n)
    X <- matrix(0, 3, n)                             # rows: mu0*, mu1*, delta
    X[3, ] <- ybar
    X[1, ] <- max(mean(y[S == 0]) - ybar, 1e-3); X[2, ] <- min(mean(y[S == 1]) - ybar, -1e-3)
  } else {
    S <- start$S; theta <- start$theta; s2w <- start$s2w; s2eta <- start$s2eta
    q <- start$q; p <- start$p; h <- start$h; X <- start$X
  }
  lam <- if (!is.null(start$lam)) start$lam else rep(1, n)
  nu_s <- if (is.null(prior$nu_s)) 2 else prior$nu_s; h_s <- if (is.null(prior$h_s)) 1 else prior$h_s       # constant-variance prior sigma2 ~ IG(nu_s/2, h_s/2)
  nswap <- 0
  x0 <- c(init_mu, delta0); P0 <- diag(c(mu_sd^2, mu_sd^2, delta_sd^2))
  keep <- ndraw - burn
  out <- list(Sprob = numeric(n), m0 = numeric(n), m1 = numeric(n), h = numeric(n),
              delta = numeric(n), mu0 = numeric(n), mu1 = numeric(n),     # posterior-mean paths of delta_t, mu0*_t, mu1*_t
              theta = matrix(NA, keep, 2), s2w = matrix(NA, keep, 2), s2eta = numeric(keep),
              P = matrix(NA, keep, 2))
  Hs <- function(s) c(1 - s, s, 1)

  for (it in 1:ndraw) {
    pi1 <- (1 - q) / ((1 - q) + (1 - p)); pi0 <- 1 - pi1
    # ---- Step 1: mu*_t (and delta_t) | S, params, y  (Carter-Kohn) -------
    sw  <- c(0, as.numeric(S[-1] != S[-n]))
    d10 <- c(0, as.numeric(S[-n] == 1 & S[-1] == 0))
    d01 <- c(0, as.numeric(S[-n] == 0 & S[-1] == 1))
    Fm <- function(t) {
      s <- sw[t]
      rbind(c(1 + theta[1] * pi0 * s, theta[1] * pi1 * s, 0),
            c(theta[2] * pi0 * s, 1 + theta[2] * pi1 * s, 0),
            c(0, 0, 1))
    }
    Om <- function(t) diag(c(d10[t] * s2w[1], d01[t] * s2w[2], sig2_eps) + kappa)
    if (use_cpp && !paper_rules) X <- eokim_ck_cpp(y, as.integer(S), h - log(lam), theta, pi0, pi1, s2w, sig2_eps, kappa, x0, diag(P0), signs) else {
    xf <- matrix(0, 3, n); Pf <- array(0, c(3, 3, n))
    for (t in 1:n) {
      if (t == 1) { xp <- x0; Pp <- P0 } else { F <- Fm(t); xp <- F %*% xf[, t - 1]; Pp <- F %*% Pf[, , t - 1] %*% t(F) + Om(t) }
      H <- Hs(S[t]); PH <- Pp %*% H; f <- as.numeric(t(H) %*% PH) + exp(h[t]) / lam[t]
      K <- PH / f
      xf[, t] <- xp + K * (y[t] - sum(H * xp))
      Pf[, , t] <- Pp - K %*% t(PH)
      Pf[, , t] <- (Pf[, , t] + t(Pf[, , t])) / 2
    }
    drawx <- function(m, V) {
      V <- V + diag(1e-10, 3)
      Lc <- tryCatch(t(chol(V)), error = function(e) NULL)
      if (is.null(Lc)) return(m)
      if (!signs) return(as.numeric(m + Lc %*% rnorm(3)))
      for (k in 1:30) { x <- as.numeric(m + Lc %*% rnorm(3)); if (x[1] > 0 && x[2] < 0) return(x) }
      c(max(x[1], 1e-6), min(x[2], -1e-6), x[3])
    }
    X[, n] <- drawx(xf[, n], Pf[, , n])
    for (t in (n - 1):1) {
      F <- Fm(t + 1); Pt <- Pf[, , t]
      A <- F %*% Pt %*% t(F) + Om(t + 1)
      G <- Pt %*% t(F) %*% solve(A)
      m <- xf[, t] + G %*% (X[, t + 1] - F %*% xf[, t])
      V <- Pt - G %*% F %*% Pt; V <- (V + t(V)) / 2
      if (paper_rules) {
        m <- as.numeric(m); s0 <- S[t]; s1 <- S[t + 1]
        cond_draw <- function(fixed) {              # fix one mu component at its conditional mean, draw the other two given it
          free <- setdiff(1:3, fixed); Vf <- V[free, free] - V[free, fixed, drop = FALSE] %*% V[fixed, free, drop = FALSE] / max(V[fixed, fixed], 1e-12)
          Vf <- (Vf + t(Vf)) / 2 + diag(1e-12, 2); Lc <- t(chol(Vf)); x <- m
          for (k in 1:30) { x[free] <- m[free] + as.numeric(Lc %*% rnorm(2)); if (!signs || (x[1] > 0 && x[2] < 0)) break }
          if (signs) { x[1] <- max(x[1], 1e-6); x[2] <- min(x[2], -1e-6) }; x }
        X[, t] <- if (s0 == 0 && s1 == 1) cond_draw(1) else if (s0 == 1 && s1 == 0) cond_draw(2) else {
          x <- m; x[3] <- m[3] + sqrt(max(V[3, 3], 0)) * rnorm(1); if (signs) { x[1] <- max(x[1], 1e-6); x[2] <- min(x[2], -1e-6) }; x }
      } else X[, t] <- drawx(as.numeric(m), V)
    }
    }
    mu0 <- X[1, ]; mu1 <- X[2, ]; dl <- X[3, ]

    # ---- Step 2: S_t | S_{-t}, mu*, h  (single move, Albert-Chib style) ---
    Tm <- matrix(c(q, 1 - q, 1 - p, p), 2, 2, byrow = TRUE)
    sd_t <- exp((h - log(lam)) / 2)
    l0 <- dnorm(y, dl + mu0, sd_t); l1 <- dnorm(y, dl + mu1, sd_t)
    if (update_S && use_cpp) S <- eokim_S_cpp(y, as.integer(S), dl, mu0, mu1, h - log(lam), q, p, pi0, pi1, if (is.null(init_p1)) -1 else init_p1) else if (update_S) for (t in 1:n) {
      pr <- if (t == 1) { if (is.null(init_p1)) c(pi0, pi1) else c(1 - init_p1, init_p1) } else Tm[S[t - 1] + 1, ]
      if (t < n) pr <- pr * Tm[, S[t + 1] + 1]
      w0 <- pr[1] * l0[t]; w1 <- pr[2] * l1[t]
      tot <- w0 + w1
      S[t] <- if (!is.finite(tot) || tot <= 0) S[t] else as.integer(runif(1) < w1 / tot)
    }
    if (guard && (all(S == 0) || all(S == 1))) S[which.min(y)] <- 1L - S[1]
    if (t_errors) { e <- y - dl - ifelse(S == 1, mu1, mu0); lam <- rgamma(n, (nu + 1) / 2, rate = (nu + e^2 / exp(h)) / 2) }        # lam_t | S, X, h

    # ---- Step 3: theta_j, sigma2_omega_j from episode shocks -------------
    sw  <- c(0, as.numeric(S[-1] != S[-n]))
    ti <- list(which(c(FALSE, S[-n] == 1 & S[-1] == 0)),    # onsets of booms  (j = 0)
               which(c(FALSE, S[-n] == 0 & S[-1] == 1)))    # onsets of recessions (j = 1)
    pi1 <- (1 - q) / ((1 - q) + (1 - p)); pi0 <- 1 - pi1
    for (j in 1:2) {
      tt <- ti[[j]]; Nj <- length(tt)
      xj <- if (Nj) pi0 * mu0[tt - 1] + pi1 * mu1[tt - 1] else numeric(0)
      dj <- if (Nj) X[j, tt] - X[j, tt - 1] else numeric(0)
      if (update_theta) for (rep in 1:50) {
        prec <- 1 / prior$theta_v + sum(xj^2) / s2w[j]
        mbar <- (prior$theta_m / prior$theta_v + sum(xj * dj) / s2w[j]) / prec
        th <- rnorm(1, mbar, sqrt(1 / prec))
        thn <- theta; thn[j] <- th
        if (abs(1 + thn[1] * pi0 + thn[2] * pi1) < 1) { theta <- thn; break }
      }
      om <- dj - theta[j] * xj
      s2w[j] <- 1 / rgamma(1, (prior$nu_w + Nj) / 2, rate = (prior$h_w + sum(om^2)) / 2)
    }

    # ---- Step 4: stochastic volatility h_t (KSC mixture) ------------------
    if (const_var) {
      e <- sqrt(lam) * (y - dl - ifelse(S == 1, mu1, mu0)); s2e <- 1 / rgamma(1, (nu_s + n) / 2, rate = (h_s + sum(e^2)) / 2); h <- rep(log(s2e), n)
    } else if (update_h) {
      e <- sqrt(lam) * (y - dl - ifelse(S == 1, mu1, mu0))
      h <- if (sv_single) sv_single_move(e, h, s2eta, h0m, h0v) else sv_draw(e, h, if (sv_sd_mode) s2eta^2 else s2eta, h0m, h0v, mix = sv_mix, cc = sv_c, bug = sv_bug)
      if (sv_prior_sd) { ssq <- sum(diff(h)^2); grd <- exp(seq(log(0.005), log(3), length.out = 600)); a_ <- prior$nu_eta / 2; b_ <- prior$h_eta / 2
        lp <- -(n - 1) * log(grd) - ssq / (2 * grd^2) - (a_ + 1) * log(grd) - b_ / grd + log(grd)          # last term: grid is uniform in log(sigma)
        s2eta <- sample(grd, 1, prob = exp(lp - max(lp)))^2
      } else s2eta <- 1 / rgamma(1, (prior$nu_eta + n - 1) / 2, rate = (prior$h_eta + sum(diff(switch(s2eta_incr, log = h, level = exp(h), sd = exp(h / 2)))^2)) / 2)
    }

    # ---- Step 5: transition probabilities -------------------------------
    N00 <- sum(S[-n] == 0 & S[-1] == 0); N01 <- sum(S[-n] == 0 & S[-1] == 1)
    N10 <- sum(S[-n] == 1 & S[-1] == 0); N11 <- sum(S[-n] == 1 & S[-1] == 1)
    if (update_qp) for (rep in 1:500) {
      qn <- rbeta(1, prior$B00 + N00, prior$B01 + N01)
      pn <- rbeta(1, prior$B11 + N11, prior$B10 + N10)
      pn1 <- (1 - qn) / ((1 - qn) + (1 - pn))
      if (!qp_constraint || abs(1 + theta[1] * (1 - pn1) + theta[2] * pn1) < 1) { q <- qn; p <- pn; break }
    }

    if (it > burn) {
      k <- it - burn
      swp <- !signs && mean(mu0 - mu1) < 0                    # no sign restrictions: store the draw in the identified labelling (booms = the higher-mean regime)
      if (swp) nswap <- nswap + 1
      Sk <- if (swp) 1L - S else S; m0k <- if (swp) mu1 else mu0; m1k <- if (swp) mu0 else mu1
      out$Sprob <- out$Sprob + Sk / keep
      out$m0 <- out$m0 + (dl + m0k) / keep; out$m1 <- out$m1 + (dl + m1k) / keep
      out$h  <- out$h + h / keep
      out$delta <- out$delta + dl / keep; out$mu0 <- out$mu0 + m0k / keep; out$mu1 <- out$mu1 + m1k / keep
      out$theta[k, ] <- if (swp) rev(theta) else theta; out$s2w[k, ] <- if (swp) rev(s2w) else s2w; out$s2eta[k] <- s2eta; out$P[k, ] <- if (swp) c(p, q) else c(q, p)
    }
    if (verbose && it %% 250 == 0) message("Eo-Kim draw ", it, "/", ndraw)
  }
  out$swap_share <- nswap / keep
  out$state <- list(S = S, theta = theta, s2w = s2w, s2eta = s2eta, q = q, p = p, h = h, X = X, lam = lam)
  out
}
