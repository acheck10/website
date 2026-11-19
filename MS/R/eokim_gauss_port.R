# Line-by-line R port of Eo & Kim's own GAUSS code for Model II (eo_kim_replication_code/m2.g + procedures.g), with switches so each feature that differs from my own
# implementation (eokim_gibbs) can be turned off one at a time, and a one-sweep function for Geweke tests.
#
# What their code does (defaults of gp_defaults()):
#   * priors actually used: sigma2_eta ~ IG(5, .05) (SIG2H_PRIOR = 10|10*0.01), NOT the printed IG(4,1); variance sigma2_e,t = B * exp(h_t) with h_1 FIXED at 0, h a random walk,
#     B ~ IG(5,5); sigma2_eps ~ IG(5000,2) (updated; essentially .0004); sigma2_omega ~ IG(5,1); theta ~ N(-.1,.25); q ~ Beta(9,1), p ~ Beta(8,2);
#     initial state (mu0*, mu1*, delta) ~ N((1,-1,.9), I) at t = 0 (the first period carries a switch)
#   * S drawn JOINTLY (Hamilton filter + backward sampling) given the mu* path; S_1 = S_2 = 0 forced; a draw with mean(S) < .2 or > .8 is redrawn
#   * chain started at the NBER partition (S = NBER, delta = mean(y), theta = -.5, sig = .3, q = .95, p = .8)
#   * mu* path: Carter-Kohn with the paper's backward rules (i)-(iii) (mu* set to conditional means in non-switch periods) and sign restrictions mu0* > 0, mu1* < 0 (<= 1000 tries)
#   * single-move accept/reject draw of the log-volatility path
# Options (opts list): pin_h1, s_draw ("joint"|"single"), s_rule (c(lo,hi) or NULL), force_first (0..), start ("nber"|"quantile"), plus all prior hyperparameters (see gp_defaults).
.gam_ig <- function(v0, d0, ee) (d0 + sum(ee^2)) / 2 / rgamma(1, (v0 + length(ee)) / 2)     # GAUSS gen_sigma2 / IG: IG((v0+n)/2, (d0+SS)/2)
.d_switch <- function(s) { n <- length(s); d10 <- d01 <- numeric(n); for (i in 2:n) { if (s[i - 1] == 1 && s[i] == 0) d10[i] <- 1; if (s[i - 1] == 0 && s[i] == 1) d01[i] <- 1 }; list(d10 = d10, d01 = d01) }

gp_defaults <- function() list(
  Pcap_ = rbind(c(9, 1), c(2, 8)), th_m = -0.1, th_v = 0.25, v1 = 10, v2 = 2, om1 = 10000, om2 = 4,
  sig2h_prior = c(10, 0.1), sig2b_prior = c(10, 10), h0 = 0, init_state = c(1, -1, 0.9),
  pin_h1 = TRUE, h1_prior = c(NA, 5),            # pin_h1 = FALSE: B = 1 and h_1 free with prior N(h1_prior[1] (NA: log var(y)), h1_prior[2])
  s_draw = "joint", s_rule = c(0.2, 0.8), force_first = 2, start = "nber", q0 = 0.95, p0 = 0.8,
  mu_rule = "paper",                 # "paper": non-switch periods set mu0*, mu1* to conditional means (rules i-iii); "carry": they are carried back unchanged within the episode (as the model implies)
  max_redraw = 2000,                 # guard for the S-rule rejection loop (the GAUSS code has none)
  sv_max_tries = 100)                # guard for the single-move accept/reject loop: after this many rejections do one Metropolis step with the prior proposal instead (the KSC proposal
                                     # essentially never accepts for a huge outlier such as 2020, and the GAUSS loop would not terminate)

gp_init <- function(y0, nber, o) {
  n <- length(y0); s1 <- if (o$start == "nber") as.integer(nber) else as.integer(y0 < quantile(y0, 0.15))
  if (o$force_first > 0) s1[1:o$force_first] <- 0L
  dd <- .d_switch(s1); D10 <- dd$d10; D01 <- dd$d01; if (s1[1] == 0) D10[1] <- 1 else D01[1] <- 1
  Pcap <- rbind(c(o$q0, 1 - o$q0), c(1 - o$p0, o$p0))
  list(s1 = s1, D10 = D10, D01 = D01, delta = rep(mean(y0), n), mu = cbind(rep(0.3, n), rep(-0.3, n)), sig_d = 0.001, theta0 = -0.5, theta1 = -0.5, sig0 = 0.3, sig1 = 0.3, Pcap = Pcap,
       sig2h = 0.1^2, H = rep(o$h0, n), sig2b = 1, sigmat = rep(0.8, n)) }

gp_un <- function(P) { un1 <- (1 - P[1, 1]) / (2 - P[1, 1] - P[2, 2]); c(un0 = 1 - un1, un1 = un1) }

gp_sweep <- function(y0, st, o) {
  n <- length(y0); un <- gp_un(st$Pcap); un0 <- un[["un0"]]; un1 <- un[["un1"]]
  ch <- function(P) tryCatch(t(chol(P + diag(1e-12, nrow(P)))), error = function(e) diag(sqrt(pmax(diag(P), 1e-12)), nrow(P)))
  # ---- 1. mu* and delta paths (Carter-Kohn, backward rules (i)-(iii), sign restrictions)
  Fm <- function(t) { d <- st$D10[t] + st$D01[t]; rbind(c(1 + st$theta0 * un0 * d, st$theta0 * un1 * d, 0), c(st$theta1 * un0 * d, 1 + st$theta1 * un1 * d, 0), c(0, 0, 1)) }
  Att <- matrix(0, n, 3); Ptt <- array(0, c(n, 3, 3)); A1 <- o$init_state; P1 <- diag(3)
  for (t in 1:n) {
    F <- Fm(t); H <- c(1 - st$s1[t], st$s1[t], 1); Q <- diag(c(st$D10[t] * st$sig0, st$D01[t] * st$sig1, st$sig_d)); R <- st$sigmat[t]
    a1 <- F %*% A1; p1 <- F %*% P1 %*% t(F) + Q; p1 <- 0.5 * (p1 + t(p1)); e <- y0[t] - sum(H * a1); f <- as.numeric(t(H) %*% p1 %*% H) + R
    A1 <- as.numeric(a1 + p1 %*% H / f * e); P1 <- p1 - (p1 %*% H) %*% t(H) %*% p1 / f; P1 <- 0.5 * (P1 + t(P1)); Att[t, ] <- A1; Ptt[t, , ] <- P1 }
  mu3 <- matrix(0, n, 3); tmp <- NULL
  for (jj in 1:1001) { tmp <- Att[n, ] + ch(Ptt[n, , ]) %*% rnorm(3); if (tmp[1] > 0 && tmp[2] < 0) break }
  ms <- as.numeric(tmp); mu3[n, ] <- ms
  for (t in (n - 1):1) {
    Q <- diag(c(st$D10[t + 1] * st$sig0, st$D01[t + 1] * st$sig1, st$sig_d)); F <- Fm(t + 1); At_ <- Att[t, ]; Pt_ <- Ptt[t, , ]
    if (st$D10[t + 1] == 1) {
      pm <- c(2, 1, 3); ms <- ms[pm]; Q <- Q[pm, pm]; F <- F[pm, pm]; At_ <- At_[pm]; Pt_ <- Pt_[pm, pm]
      F23 <- F[2:3, , drop = FALSE]; P1t <- F23 %*% Pt_ %*% t(F23) + Q[2:3, 2:3]; P1t <- 0.5 * (P1t + t(P1t)); G <- Pt_ %*% t(F23) %*% solve(P1t)
      At <- as.numeric(At_ + G %*% (ms[2:3] - F23 %*% At_)); Pt <- Pt_ - G %*% F23 %*% Pt_; Pt <- 0.5 * (Pt + t(Pt)); ms <- At
      for (jj in 1:1001) { tmp <- At[2:3] + ch(Pt[2:3, 2:3]) %*% rnorm(2); if (tmp[1] > 0) { ms[2:3] <- tmp; break } }
      ms <- ms[c(2, 1, 3)]
    } else if (st$D01[t + 1] == 1) {
      F23 <- F[2:3, , drop = FALSE]; P1t <- F23 %*% Pt_ %*% t(F23) + Q[2:3, 2:3]; P1t <- 0.5 * (P1t + t(P1t)); G <- Pt_ %*% t(F23) %*% solve(P1t)
      At <- as.numeric(At_ + G %*% (ms[2:3] - F23 %*% At_)); Pt <- Pt_ - G %*% F23 %*% Pt_; Pt <- 0.5 * (Pt + t(Pt)); ms <- At
      for (jj in 1:1001) { tmp <- At[2:3] + ch(Pt[2:3, 2:3]) %*% rnorm(2); if (tmp[1] < 0) { ms[2:3] <- tmp; break } }
    } else {
      F3 <- F[3, , drop = FALSE]; P1t <- as.numeric(F3 %*% Pt_ %*% t(F3)) + st$sig_d; G <- Pt_ %*% t(F3) / P1t
      At <- as.numeric(At_ + G * (ms[3] - sum(F3 * At_))); Pt <- Pt_ - G %*% F3 %*% Pt_; Pt <- 0.5 * (Pt + t(Pt)); keep12 <- ms[1:2]; ms <- At; ms[3] <- At[3] + sqrt(max(Pt[3, 3], 0)) * rnorm(1)
      if (o$mu_rule == "carry") ms[1:2] <- keep12 }
    mu3[t, ] <- ms }
  mu <- mu3[, 1:2]; delta <- mu3[, 3]; st$mu <- mu; st$delta <- delta
  # ---- 2. sigma2_eps, 3. volatility path, sigma2_eta, scale B
  st$sig_d <- .gam_ig(o$om1, o$om2, diff(delta))
  must <- mu[, 1] * (1 - st$s1) + mu[, 2] * st$s1; yy <- y0 - delta; E <- y0 - must - delta
  H <- st$H
  sv_step <- function(dta, h, s2h, free_first) {
    first <- if (free_first) 1 else 2; if (!free_first) h[1] <- o$h0
    for (ii in first:n) {
      if (ii == 1) { h1m <- if (is.na(o$h1_prior[1])) log(var(y0 - delta)) else o$h1_prior[1]; pr <- 1 / o$h1_prior[2] + 1 / s2h; hs <- (h1m / o$h1_prior[2] + h[2] / s2h) / pr; s2 <- 1 / pr
      } else if (ii == n) { hs <- h[ii - 1]; s2 <- s2h } else { hs <- 0.5 * (h[ii - 1] + h[ii + 1]); s2 <- 0.5 * s2h }
      mu_ <- hs + 0.5 * s2 * ((dta[ii]^2) * exp(-hs) - 1)
      ntry <- 0
      repeat { ntry <- ntry + 1; ht <- mu_ + sqrt(s2) * rnorm(1); acnu <- exp(-0.5 * dta[ii]^2 * exp(-ht)); acdn <- exp(-0.5 * dta[ii]^2 * (exp(-hs) * (1 + hs) - ht * exp(-hs))); acc <- acnu / acdn
        if (is.na(acc)) acc <- 0
        if (runif(1) < acc) { h[ii] <- ht; break }
        if (ntry >= o$sv_max_tries) { hp <- hs + sqrt(s2) * rnorm(1); llp <- -0.5 * hp - 0.5 * dta[ii]^2 * exp(-hp); llc <- -0.5 * h[ii] - 0.5 * dta[ii]^2 * exp(-h[ii]); if (log(runif(1)) < llp - llc) h[ii] <- hp; break } } }
    h }
  if (o$pin_h1) { H <- sv_step(E / sqrt(st$sig2b), H, st$sig2h, FALSE) } else { st$sig2b <- 1; H <- sv_step(E, H, st$sig2h, TRUE) }
  st$sig2h <- .gam_ig(o$sig2h_prior[1], o$sig2h_prior[2], diff(H))
  if (o$pin_h1) st$sig2b <- .gam_ig(o$sig2b_prior[1], o$sig2b_prior[2], E / sqrt(exp(H)))
  st$H <- H; st$sigmat <- st$sig2b * exp(H)
  # ---- 4. regime path S
  q <- st$Pcap[1, 1]; p <- st$Pcap[2, 2]; sg <- st$sigmat
  f0 <- function(t) exp(-0.5 * (yy[t] - mu[t, 1])^2 / sg[t]) / sqrt(2 * pi * sg[t]); f1 <- function(t) exp(-0.5 * (yy[t] - mu[t, 2])^2 / sg[t]) / sqrt(2 * pi * sg[t])
  ff <- o$force_first
  if (o$s_draw == "joint") {
    ntry <- 0
    repeat {
      ntry <- ntry + 1; flt <- matrix(0, 2, n); l0 <- (1 - p) / (2 - p - q); l1 <- (1 - q) / (2 - p - q)
      for (t in 1:n) { t0 <- q * l0 + (1 - p) * l1; t1 <- (1 - q) * l0 + p * l1; la0 <- -0.5 * (yy[t] - mu[t, 1])^2 / sg[t]; la1 <- -0.5 * (yy[t] - mu[t, 2])^2 / sg[t]; mx <- max(la0, la1); a0 <- exp(la0 - mx); a1 <- exp(la1 - mx); fy <- a0 * t0 + a1 * t1; l0 <- a0 * t0 / fy; l1 <- a1 * t1 / fy; flt[, t] <- c(l0, l1) }
      s <- integer(n); s[n] <- as.integer(runif(1) >= flt[1, n] / (flt[1, n] + flt[2, n]))
      for (j in (n - 1):(ff + 1)) { if (s[j + 1] == 0) { p0 <- q * flt[1, j]; p1 <- (1 - p) * flt[2, j] } else { p0 <- (1 - q) * flt[1, j]; p1 <- p * flt[2, j] }
        s[j] <- as.integer(runif(1) >= p0 / (p0 + p1)) }
      if (is.null(o$s_rule) || (mean(s) >= o$s_rule[1] && mean(s) <= o$s_rule[2]) || ntry >= o$max_redraw) break }
  } else {                                                                # single-move (Albert-Chib) draw given the mu* path
    s <- st$s1; pi1 <- un1; Tm <- rbind(c(q, 1 - q), c(1 - p, p))
    for (t in (ff + 1):n) { pr <- if (t == 1) c(1 - pi1, pi1) else Tm[s[t - 1] + 1, ]; if (t < n) pr <- pr * Tm[, s[t + 1] + 1]
      w0 <- pr[1] * f0(t); w1 <- pr[2] * f1(t); s[t] <- if (w0 + w1 > 0) as.integer(runif(1) < w1 / (w0 + w1)) else s[t] }
    if (!is.null(o$s_rule) && (mean(s) < o$s_rule[1] || mean(s) > o$s_rule[2])) s <- st$s1 }
  if (ff > 0) s[1:ff] <- 0L
  dd <- .d_switch(s); st$D10 <- dd$d10; st$D01 <- dd$d01; if (s[1] == 0) st$D10[1] <- 1 else st$D01[1] <- 1; st$s1 <- s
  # ---- 5. transition probabilities
  cnt <- matrix(0, 2, 2); for (j in 2:n) cnt[s[j - 1] + 1, s[j] + 1] <- cnt[s[j - 1] + 1, s[j] + 1] + 1; sw <- o$Pcap_ + cnt
  pj <- c(rbeta(1, sw[1, 1], sw[1, 2]), rbeta(1, sw[2, 1], sw[2, 2])); st$Pcap <- rbind(c(pj[1], 1 - pj[1]), c(pj[2], 1 - pj[2]))
  un <- gp_un(st$Pcap); un0 <- un[["un0"]]; un1 <- un[["un1"]]
  # ---- 6. theta_j, sigma2_omega_j (re-drawn until the stationarity bound holds)
  m0p <- mu[st$D10 == 1, 1]; m0sp <- mu[st$D01 == 1, 1]; m1p <- mu[st$D01 == 1, 2]; m1sp <- mu[st$D10 == 1, 2]
  gb <- function(b0, B0, sg_, yv, xv) { B1 <- 1 / (1 / B0 + sum(xv^2) / sg_); b1 <- B1 * (b0 / B0 + sum(xv * yv) / sg_); rnorm(1, b1, sqrt(B1)) }
  repeat {
    L0 <- max(0, min(length(m0p) - 1, length(m0sp), length(m1p))); i0 <- seq_len(L0)                  # = their k / k-1 on real data; also safe when there are very few episodes
    y0v <- m0p[i0 + 1] - m0sp[i0]; x0v <- m0sp[i0] * un0 + m1p[i0] * un1
    th0 <- gb(o$th_m, o$th_v, st$sig0, y0v, x0v); sg0 <- .gam_ig(o$v1, o$v2, y0v - x0v * th0)
    L1 <- max(0, min(length(m1p), length(m1sp), length(m0p))); i1 <- seq_len(L1); y1v <- m1p[i1] - m1sp[i1]; x1v <- m0p[i1] * un0 + m1sp[i1] * un1
    th1 <- gb(o$th_m, o$th_v, st$sig1, y1v, x1v); sg1 <- .gam_ig(o$v1, o$v2, y1v - th1 * x1v)
    st$theta0 <- th0; st$sig0 <- sg0; st$theta1 <- th1; st$sig1 <- sg1
    if (abs(1 + th0 * un0 + th1 * un1) < 1) break }
  st }

eokim_gauss_port <- function(y0, nber, ndraw = 12000, burn = 4000, seed = 1, verbose = FALSE, opts = list()) {
  set.seed(seed); n <- length(y0); o <- modifyList(gp_defaults(), opts); st <- gp_init(y0, nber, o); keep <- ndraw - burn
  out <- list(s1 = numeric(n), sv = numeric(n), h = numeric(n), must_y = numeric(n), delta = numeric(n), mu0 = numeric(n), mu1 = numeric(n),
              par = matrix(NA_real_, keep, 8, dimnames = list(NULL, c("theta0", "theta1", "sig0", "sig1", "sig_d", "s2eta", "q", "p"))), s2b = numeric(keep))
  for (it in 1:ndraw) {
    st <- gp_sweep(y0, st, o)
    if (it > burn) { k <- it - burn; must <- st$mu[, 1] * (1 - st$s1) + st$mu[, 2] * st$s1
      out$par[k, ] <- c(st$theta0, st$theta1, st$sig0, st$sig1, st$sig_d, st$sig2h, st$Pcap[1, 1], st$Pcap[2, 2]); out$s2b[k] <- st$sig2b
      out$s1 <- out$s1 + st$s1 / keep; out$sv <- out$sv + st$sigmat / keep; out$h <- out$h + st$H / keep; out$must_y <- out$must_y + (must + st$delta) / keep; out$delta <- out$delta + st$delta / keep
      out$mu0 <- out$mu0 + (st$mu[, 1] + st$delta) / keep; out$mu1 <- out$mu1 + (st$mu[, 2] + st$delta) / keep }
    if (verbose && it %% 500 == 0) message("iteration ", it, "/", ndraw) }
  out$state <- st; out
}

# ---- option sets for the website data -------------------------------------------------------------------------------------------------------------------------
# "quarterly": exactly their GAUSS settings (quarterly GDP growth, unit-variance scale).
# "monthly"  : monthly data (monthly common factor, simulated 700-month series). Rescaled as follows: transition priors keep the same prior MEANS in expected
#              durations (q: .967 = 30 months; p: .933 = 15 months; prior strength 30 / 30 pseudo-observations: Beta(29, 1) and Beta(28, 2), the same as trans_prior("monthly") in ms_samplers.R); the long-run-mean random-walk variance and the
#              volatility innovation variance are divided by 3 (quarterly -> monthly: same variation per year); the scale B and the initial regime means are
#              set for the factor's scale (variance ~ .4, boom/recession deviations ~ +.3 / -1); the 20-80% recession-share rule becomes 3-50% for the simulations and 14-80% for the real monthly factor (their 20% lower bound is 1.2x
#              the NBER quarterly share of 17%; NBER recessions are 12% of months, hence 14%); chain start values q0 = .98, p0 = .9.
gp_opts <- function(kind = c("quarterly", "monthly"), y = NULL, start = "nber", s_rule = c(0.03, 0.5)) {
  kind <- match.arg(kind)
  if (kind == "quarterly") return(list(start = start))
  list(start = start, Pcap_ = rbind(c(29, 1), c(2, 28)), om2 = 1.33, sig2h_prior = c(10, 0.0333), sig2b_prior = c(10, 3), init_state = c(0.3, -1, if (is.null(y)) 0 else mean(y)),
       s_rule = s_rule, q0 = 0.98, p0 = 0.9)
}
