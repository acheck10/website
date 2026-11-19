if (!exists(".src_cpp")) source("R/rcpp_cache.R")      # persistent, self-repairing compile cache for the C++ kernels
# Exact Gibbs / Metropolis-within-Gibbs sampler for the Eo & Kim (2013) model, with the
# regime-specific mean path mu*_t integrated out where the paper's sampler conditions on it.
#
# Model: see ms_models_extra.R (extended model with random-walk long-run growth and stochastic volatility).
# Differences from eokim_gibbs() (the replication of the paper's algorithm):
#   * S is drawn single-move from p(S_t | S_-t, y, params) with the mu* path INTEGRATED OUT (Kalman likelihood,
#     C++ in eokim_kf.cpp). The paper's step conditions on a mu* path that itself depends on S.
#   * theta_j, sigma2_omega_j and (q, p) are updated by Metropolis-Hastings on the same collapsed likelihood
#     (random-walk proposals for theta and log s2w; an independence proposal from the transition-count Beta
#     posterior for (q, p), which also carries the Pr(S_1 | q, p) factor and the stationarity restriction).
#   * mu*, delta path | S, params, h, y by Carter-Kohn; h, sigma2_eta | path by the 10-component SV draw.
#   * NO hard sign restrictions (mu0* > 0, mu1* < 0): the truncation probability of a Gaussian path has no
#     closed form, so an exact sampler cannot impose them. Identification relies on the prior x_1 ~ N((.5,-.5,.), I).
#
# Needs: ms_models_extra.R (eokim_prior, sv_draw, .omori10), Rcpp.

suppressPackageStartupMessages(library(Rcpp))
if (!exists("kf_loglik")) .src_cpp("R/eokim_kf.cpp")
if (!exists("eokim_ck_cpp")) .src_cpp("R/hamilton.cpp")

# Carter-Kohn draw of x_t = (mu0*, mu1*, delta)_t given S, parameters, h and y.
# (State noise is singular at non-switch dates; a tiny ridge kappa regularises the backward step.)
ck_draw_X_R <- function(y, S, theta, s2w, pi0, pi1, h, x0, P0, sig2_eps, kappa = 1e-9) {      # reference R implementation (ck_draw_X below calls the C++ port, same algorithm / random-number order)
  n <- length(y)
  sw  <- c(0, as.numeric(S[-1] != S[-n]))
  d10 <- c(0, as.numeric(S[-n] == 1 & S[-1] == 0))
  d01 <- c(0, as.numeric(S[-n] == 0 & S[-1] == 1))
  Fm <- function(t) { s <- sw[t]
    rbind(c(1 + theta[1] * pi0 * s, theta[1] * pi1 * s, 0), c(theta[2] * pi0 * s, 1 + theta[2] * pi1 * s, 0), c(0, 0, 1)) }
  Om <- function(t) diag(c(d10[t] * s2w[1], d01[t] * s2w[2], sig2_eps) + kappa)
  xf <- matrix(0, 3, n); Pf <- array(0, c(3, 3, n))
  for (t in 1:n) {
    if (t == 1) { xp <- x0; Pp <- P0 } else { F <- Fm(t); xp <- F %*% xf[, t - 1]; Pp <- F %*% Pf[, , t - 1] %*% t(F) + Om(t) }
    H <- c(1 - S[t], S[t], 1); PH <- Pp %*% H; f <- as.numeric(t(H) %*% PH) + exp(h[t]); K <- PH / f
    xf[, t] <- xp + K * (y[t] - sum(H * xp)); Pt <- Pp - K %*% t(PH); Pf[, , t] <- (Pt + t(Pt)) / 2
  }
  drawx <- function(m, V) {
    V <- V + diag(1e-10, 3)
    Lc <- tryCatch(t(chol(V)), error = function(e) NULL)
    if (is.null(Lc)) return(m)
    as.numeric(m + Lc %*% rnorm(3))
  }
  X <- matrix(0, 3, n)
  X[, n] <- drawx(xf[, n], Pf[, , n])
  for (t in (n - 1):1) {
    F <- Fm(t + 1); Pt <- Pf[, , t]
    A <- F %*% Pt %*% t(F) + Om(t + 1)
    G <- Pt %*% t(F) %*% solve(A)
    m <- xf[, t] + G %*% (X[, t + 1] - F %*% xf[, t])
    V <- Pt - G %*% F %*% Pt; V <- (V + t(V)) / 2
    X[, t] <- drawx(as.numeric(m), V)
  }
  X
}

ck_draw_X <- function(y, S, theta, s2w, pi0, pi1, h, x0, P0, sig2_eps, kappa = 1e-9)
  eokim_ck_cpp(y, as.integer(S), h, theta, pi0, pi1, s2w, sig2_eps, kappa, x0, diag(P0), FALSE)

eokim_exact_gibbs <- function(y, ndraw = 3000, burn = 1000, prior = eokim_prior,
                              sig2_eps = 0.0004, init_mu = c(0.5, -0.5), kappa = 1e-9,
                              seed = NULL, verbose = TRUE, start = NULL,
                              delta0 = NULL, h0m = NULL, h0v = 5, init_p1 = NULL,
                              sv_mix = .omori10, sv_c = 1e-8,
                              step_theta = 0.35, step_lw = 0.8, n_mh = 3, adapt = TRUE, init_swap = FALSE, delta_sd = 1, mu_sd = 1, restrict = NULL, max_x_tries = 100000, n_block = 0L, block_shift = 6L, block_win = 12L) {
  # n_block > 0: after each single-move sweep, n_block Metropolis block moves on the regime path (boundary shifts of up to block_shift months, window flips of up to block_win months; sm_block_S in eokim_kf.cpp), same collapsed target.
  # restrict = c(c0, c1, c2): impose  c' mean_t(x_t) > 0  EXACTLY on the time average of the path x_t = (mu0*, mu1*, delta)  (joint truncation of the prior: the prior is conditioned on the event).
  # Presets: c(1, -1, 0)  mean(mu0*) > mean(mu1*) (labels ordered: no label switching, no relabelling needed);  c(0, -1, 0)  mean(mu1*) < 0 (the analogue of the low-regime-mean < 0 restriction).
  # The collapsed posterior is multiplied by Pr(c' xbar > 0 | y, S, params, h), a Gaussian probability available in closed form from the Kalman filter with the running sum appended to the state (eokim_kf.cpp);
  # the path is drawn by Carter-Kohn with rejection until the restriction holds (acceptance probability = that same Pr).  The paper's restriction at EVERY date is not used: under the unrestricted posterior it holds with probability ~0.1%.
  # init_swap = TRUE starts the chain in the WRONG labelling (recession = high y); used to test the diagnostic.
  # step_theta / step_lw: random-walk MH step sizes. If adapt = TRUE they are tuned toward ~35% acceptance
  # during the burn-in only (every 50 sweeps), then frozen, so the retained draws come from a fixed kernel.
  # (In the Geweke test each call is a single sweep with burn = 0, so no adaptation happens there.)
  if (!is.null(seed)) set.seed(seed)
  n <- length(y); ybar <- mean(y)
  if (is.null(delta0)) delta0 <- ybar
  if (is.null(h0m)) h0m <- log(var(y))
  if (is.null(start)) {
    S <- if (init_swap) as.integer(y > quantile(y, 0.85)) else as.integer(y < quantile(y, 0.15))
    theta <- c(-0.1, -0.1); s2w <- c(0.25, 0.25); s2eta <- 0.05; q <- 0.9; p <- 0.8
    h <- rep(log(var(y)), n)
    X <- matrix(0, 3, n)
  } else { S <- start$S; theta <- start$theta; s2w <- start$s2w; s2eta <- start$s2eta
           q <- start$q; p <- start$p; h <- start$h; X <- start$X }
  x0 <- c(init_mu, delta0); P0 <- diag(c(mu_sd^2, mu_sd^2, delta_sd^2)); P0v <- as.numeric(t(P0))
  init1 <- if (is.null(init_p1)) -1 else init_p1
  ig_a <- prior$nu_w / 2; ig_b <- prior$h_w / 2
  pis <- function(q, p) { p1 <- (1 - q) / ((1 - q) + (1 - p)); c(1 - p1, p1) }
  stat_ok <- function(th, q, p) { pp <- pis(q, p); abs(1 + th[1] * pp[1] + th[2] * pp[2]) < 1 }
  if (!is.null(restrict)) { restrict <- as.numeric(restrict); stopifnot(length(restrict) == 3) }
  LL <- function(th, s2, qq, pp) if (is.null(restrict)) kf_loglik(y, S, h, th[1], th[2], s2[1], s2[2], qq, pp, sig2_eps, x0, P0v) else sum(kf_loglik_restr(y, S, h, th[1], th[2], s2[1], s2[2], qq, pp, sig2_eps, x0, P0v, restrict))
  keep <- ndraw - burn
  out <- list(Sprob = numeric(n), m0 = numeric(n), m1 = numeric(n), h = numeric(n),
              theta = matrix(NA, keep, 2), s2w = matrix(NA, keep, 2), swapped = logical(keep),
              s2eta = numeric(keep), P = matrix(NA, keep, 2), acc = c(th = 0, lw = 0, qp = 0))
  nacc <- c(th = 0, lw = 0, qp = 0); ntry <- c(th = 0, lw = 0, qp = 0)
  nacc_w <- nacc; ntry_w <- ntry; xtries_tot <- 0; bacc <- 0

  for (it in 1:ndraw) {
    # ---- 1. S | theta, s2w, q, p, h, y  (mu* integrated out) --------------------------------
    S <- as.integer(if (is.null(restrict)) sm_sweep_S(y, S, h, theta[1], theta[2], s2w[1], s2w[2], q, p, sig2_eps, x0, P0v, init1) else sm_sweep_S_restr(y, S, h, theta[1], theta[2], s2w[1], s2w[2], q, p, sig2_eps, x0, P0v, init1, restrict, FALSE))
    if (n_block > 0) { bm <- sm_block_S(y, S, h, theta[1], theta[2], s2w[1], s2w[2], q, p, sig2_eps, x0, P0v, init1, if (is.null(restrict)) c(0, 0, 0) else restrict, !is.null(restrict), as.integer(n_block), as.integer(block_shift), as.integer(block_win)); S <- as.integer(bm$S); bacc <- bacc + bm$acc / ndraw }
    cur <- LL(theta, s2w, q, p)
    # ---- 2. theta_j, s2w_j | S, ... (random-walk MH on the collapsed likelihood) -----------
    for (rep in 1:n_mh) {
      for (j in 1:2) {
        thn <- theta; thn[j] <- theta[j] + rnorm(1, 0, step_theta); ntry["th"] <- ntry["th"] + 1
        if (stat_ok(thn, q, p)) {
          lln <- LL(thn, s2w, q, p)
          if (log(runif(1)) < lln - cur + dnorm(thn[j], prior$theta_m, sqrt(prior$theta_v), log = TRUE) -
                              dnorm(theta[j], prior$theta_m, sqrt(prior$theta_v), log = TRUE)) {
            theta <- thn; cur <- lln; nacc["th"] <- nacc["th"] + 1 }
        }
      }
      for (j in 1:2) {
        s2n <- s2w; s2n[j] <- exp(log(s2w[j]) + rnorm(1, 0, step_lw)); ntry["lw"] <- ntry["lw"] + 1
        lln <- LL(theta, s2n, q, p)
        lpr <- function(s) -ig_a * log(s) - ig_b / s            # log IG density + log-Jacobian
        if (log(runif(1)) < lln - cur + lpr(s2n[j]) - lpr(s2w[j])) { s2w <- s2n; cur <- lln; nacc["lw"] <- nacc["lw"] + 1 }
      }
    }
    # ---- 3. (q, p) | S, ... : independence MH from the transition-count Beta posterior ---------
    N00 <- sum(S[-n] == 0 & S[-1] == 0); N01 <- sum(S[-n] == 0 & S[-1] == 1)
    N10 <- sum(S[-n] == 1 & S[-1] == 0); N11 <- sum(S[-n] == 1 & S[-1] == 1)
    qn <- rbeta(1, prior$B00 + N00, prior$B01 + N01); pn <- rbeta(1, prior$B11 + N11, prior$B10 + N10)
    ntry["qp"] <- ntry["qp"] + 1
    if (stat_ok(theta, qn, pn)) {
      lln <- LL(theta, s2w, qn, pn); lr <- lln - cur
      if (is.null(init_p1)) lr <- lr + log(pis(qn, pn)[S[1] + 1]) - log(pis(q, p)[S[1] + 1])
      if (log(runif(1)) < lr) { q <- qn; p <- pn; nacc["qp"] <- nacc["qp"] + 1 }
    }
    # ---- 4. mu*, delta path | S, params, h, y -------------------------------------------
    pp <- pis(q, p)
    X <- ck_draw_X(y, S, theta, s2w, pp[1], pp[2], h, x0, P0, sig2_eps, kappa)
    if (!is.null(restrict)) { tries <- 1; while (mean(restrict[1] * X[1, ] + restrict[2] * X[2, ] + restrict[3] * X[3, ]) <= 0) { if ((tries <- tries + 1) > max_x_tries) stop("restricted path draw: no draw satisfied the restriction in ", max_x_tries, " tries"); X <- ck_draw_X(y, S, theta, s2w, pp[1], pp[2], h, x0, P0, sig2_eps, kappa) }; xtries_tot <- xtries_tot + tries }
    # ---- 5. h, s2eta | path -------------------------------------------------------------
    e <- y - X[3, ] - ifelse(S == 1, X[2, ], X[1, ])
    h <- sv_draw(e, h, s2eta, h0m, h0v, mix = sv_mix, cc = sv_c)
    s2eta <- 1 / rgamma(1, (prior$nu_eta + n - 1) / 2, rate = (prior$h_eta + sum(diff(h)^2)) / 2)

    if (adapt && it <= burn && it %% 50 == 0) {              # burn-in step-size adaptation
      rt <- (nacc - nacc_w) / pmax(ntry - ntry_w, 1)
      step_theta <- step_theta * min(2, max(0.5, exp(2 * (rt["th"] - 0.35))))
      step_lw    <- step_lw    * min(2, max(0.5, exp(2 * (rt["lw"] - 0.35))))
      nacc_w <- nacc; ntry_w <- ntry
    }
    if (it > burn) {
      k <- it - burn
      # Label-switching diagnostic + post-hoc relabelling of the STORED draw (the chain itself is untouched,
      # so its kernel / the Geweke result are unaffected). Ordering criterion: average over t of mu0* - mu1*
      # (delta cancels). If negative, the draw sits in the swapped labelling; swap S, the mean rows,
      # (theta0,theta1), (s2w0,s2w1) and (q,p) -> (p,q).
      swapped <- if (is.null(restrict)) mean(X[1, ] - X[2, ]) < 0 else FALSE      # with a restriction the draws are in the identified labelling by construction
      out$swapped[k] <- swapped
      Sk <- if (swapped) 1L - S else S
      m0k <- X[3, ] + (if (swapped) X[2, ] else X[1, ]); m1k <- X[3, ] + (if (swapped) X[1, ] else X[2, ])
      out$Sprob <- out$Sprob + Sk / keep
      out$m0 <- out$m0 + m0k / keep; out$m1 <- out$m1 + m1k / keep; out$h <- out$h + h / keep
      out$theta[k, ] <- if (swapped) rev(theta) else theta
      out$s2w[k, ]   <- if (swapped) rev(s2w) else s2w
      out$P[k, ]     <- if (swapped) c(p, q) else c(q, p)
      out$s2eta[k] <- s2eta
    }
    if (verbose && it %% 250 == 0) message("Eo-Kim exact draw ", it, "/", ndraw)
  }
  out$xtries_mean <- xtries_tot / ndraw; out$block_acc <- bacc; out$acc <- nacc / pmax(ntry, 1); out$steps <- c(theta = step_theta, lw = step_lw)
  out$swap_share <- mean(out$swapped)    # share of retained draws found in the swapped labelling (then relabelled)
  out$state <- list(S = S, theta = theta, s2w = s2w, s2eta = s2eta, q = q, p = p, h = h, X = X)
  out
}
