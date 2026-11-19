# Multi-chain fits of every model on the DFM factor, with Gelman-style diagnostics.
# usage: Rscript scripts/run_dfm_chains.R <vintage|current> [K=4] [ndraw=4000] [burn=1000] [tag=old] [filter regex]
#   SV_VERSION env var: "old" (standard Gibbs, msmsv_gibbs/msmtsv_gibbs) or "new" (collapsed-S samplers from ms_models_sv2.R)
a <- commandArgs(trailingOnly = TRUE); vin <- a[1]
K <- if (length(a) >= 2) as.integer(a[2]) else 4L; ndraw <- if (length(a) >= 3) as.integer(a[3]) else 4000L
burn <- if (length(a) >= 4) as.integer(a[4]) else 1000L; tag <- if (length(a) >= 5) a[5] else "old"
filt <- if (length(a) >= 6) a[6] else ""
source("R/ms_samplers.R"); source("R/ms_models_extra.R"); source("R/ms_models_sv.R"); source("R/eokim_exact.R"); source("R/conv_diag.R"); source("R/priors.R")
if (tag != "old") source("R/ms_models_sv2.R")
suppressPackageStartupMessages(library(parallel))
r <- readRDS(if (grepl("^(gdp_|stress_|sim_|eokim)", vin)) sprintf("results/%s_data.rds", vin) else sprintf("results/dfm_data_%s.rds", vin));   # gdp_full / gdp_eokim: quarterly GDP growth / 4 (scripts/prep_gdp_data.R)
 y <- r$factor; S0 <- r$S; n <- length(y); keep <- ndraw - burn

make_start <- function(fam, k, a1max = NULL) {          # chain 1: sampler default; 2: NBER-informed; 3,4+: random persistent S, prior-drawn means
  if (k == 1) return(NULL)
  set.seed(1000 + k)
  if (k == 2) { S <- S0 } else { S <- integer(n); S[1] <- rbinom(1, 1, .4); for (t in 2:n) S[t] <- if (runif(1) < .93) S[t - 1] else 1L - S[t - 1] }
  if (sum(S) == 0) S[1:10] <- 1L
  a <- if (k == 2) c(mean(y[S == 0]), mean(y[S == 1])) else { a0 <- .5 + rnorm(1); a1 <- -.5 - abs(rnorm(1)); c(a0, min(a1, a0 - .1)) }
  if (!is.null(a1max)) a[2] <- min(a[2], a1max - 1e-3)
  P <- c(P00 = .93, P11 = .85); lh <- rep(log(var(y)), n)
  switch(fam, msm = list(S = S, a = a, sig2 = var(y), P = P),
              msmt = list(S = S, a = a, sig2 = var(y), lam = rep(1, n), nu = 5, P = P),
              msq = list(S = S, a = a, sig = sd(y), zs = rep(sd(y), n), P = P),
              msmsv = list(S = S, a = a, s2eta = .05, h = lh, P = P),
              msmtsv = list(S = S, a = a, s2eta = .05, h = lh, lam = rep(1, n), nu = 5, P = P), NULL)
}
sv1 <- tag == "old"
drift <- tag %in% c("drift0", "driftq"); S2D <- as.numeric(Sys.getenv("S2D", "1.4e-4"))   # tag "drift0": MSM-family models with the slow drift delta_t (fixed random-walk variance S2D) and the zero-mean split (R/ms_models_drift0.R)
if (drift) source("R/ms_models_drift0.R")
eokim_prior_tight <- modifyList(eokim_prior, list(nu_eta = 10, h_eta = 0.1, nu_w = 20, h_w = 0.9))      # replication page, exact sampler: s2eta ~ IG(5, .05) (the authors' code) and s2w ~ IG(10, .45) (mean .05: regime means move little at a switch); see investigation/EK_EXACT_RESTRICTION.md
eokim_prior_paper <- modifyList(eokim_prior, list(nu_eta = 8, h_eta = 2))      # the paper's printed prior on the volatility-innovation variance, IG(4, 1)
FREQ <- if (grepl("^(gdp_|eokim)", vin)) "quarterly" else "monthly"; TR <- trans_prior(FREQ)       # transition-probability priors follow the data frequency (R/ms_samplers.R::trans_prior); GDP is quarterly, everything else monthly
fit1 <- function(model, k) {                              # model = list(name, fam, tau, restricted)
  rs <- model$restricted; PP <- priors_for(y, FREQ); pr <- modifyList(PP$pr, if (rs) list(a1_max = 0) else list()); prs <- modifyList(PP$prs, if (rs) list(a1_max = 0) else list()); ekp <- PP$ekp; EKA <- PP$ek_args      # priors in the series' own units and calendar time (R/priors.R)
  st <- make_start(model$fam, k, if (rs) 0 else NULL); sd <- 100 * k + 7
  t0 <- proc.time()[3]
  m <- if (drift) drift0_gibbs(model$fam, y, ndraw, burn, prior = if (model$fam %in% c("msm", "msmt")) pr else prs, a1_prior = PP$drift$a1_prior, s2d = PP$drift$s2d, d0 = PP$drift$d0, dsd = PP$drift$dsd, tau = if (is.na(model$tau)) 0.3 else model$tau, seed = sd, verbose = FALSE, start = st) else switch(model$fam,
    msm = msm_gibbs(y, ndraw, burn, prior = pr, seed = sd, verbose = FALSE, start = st),
    msmt = msmt_gibbs(y, ndraw, burn, prior = pr, seed = sd, verbose = FALSE, start = st),
    msq = msq_gibbs(y, tau = model$tau, ndraw, burn, prior = pr, seed = sd, verbose = FALSE, start = st),
    msmsv = (if (sv1) msmsv_gibbs else msmsv_gibbs2)(y, ndraw, burn, prior = prs, seed = sd, verbose = FALSE, start = st),
    msmtsv = (if (sv1) msmtsv_gibbs else msmtsv_gibbs2)(y, ndraw, burn, prior = prs, seed = sd, verbose = FALSE, start = st),
    repek = eokim_gibbs(y, ndraw, burn, prior = eokim_prior_paper, seed = sd, verbose = FALSE, init_q = c(0.15, 0.10, 0.20, 0.12, 0.18, 0.25, 0.08, 0.30)[(k - 1) %% 8 + 1]),       # replication: the paper-style sampler with the paper's printed priors
    repek_ours = do.call(eokim_gibbs, c(list(y, ndraw, burn, prior = ekp, seed = sd, verbose = FALSE, signs = TRUE, init_q = c(0.15, 0.10, 0.20, 0.12, 0.18, 0.25, 0.08, 0.30)[(k - 1) %% 8 + 1]), EKA)),   # the same sampler with our standardized priors (SV, random-walk long-run mean, sign restrictions)
    repekx = eokim_exact_gibbs(y, ndraw, burn, prior = eokim_prior_tight, seed = sd, verbose = FALSE, restrict = c(0, -1, 0), n_block = as.integer(Sys.getenv("EK_NBLOCK", "0"))),      # exact sampler with the time-average restriction, tight s2w and the authors' s2eta prior (eokim_prior_tight)
    rwec = do.call(eokim_gibbs, c(list(y, ndraw, burn, prior = ekp, seed = sd, verbose = FALSE, const_var = !model$sv, t_errors = model$t, signs = model$restricted, init_q = c(0.15, 0.10, 0.20, 0.12, 0.18, 0.25, 0.08, 0.30)[(k - 1) %% 8 + 1]), modifyList(EKA, list(sig2_eps = if (model$drift) EKA$sig2_eps else 0)))),      # MSM-RW-EC: random-walk regimes with error correction (Eo-Kim-style), see below
    ek = do.call(eokim_gibbs, c(list(y, ndraw, burn, prior = ekp, seed = sd, verbose = FALSE), EKA)),
    ekx = do.call(eokim_exact_gibbs, c(list(y, ndraw, burn, prior = ekp, seed = sd, verbose = FALSE, n_block = as.integer(Sys.getenv("EK_NBLOCK", "0"))), EKA)),       # EK_NBLOCK = n: block moves for the regime path (default 0 = none)
    ekxr = do.call(eokim_exact_gibbs, c(list(y, ndraw, burn, prior = ekp, seed = sd, verbose = FALSE, restrict = c(0, -1, 0), n_block = as.integer(Sys.getenv("EK_NBLOCK", "0"))), EKA)),      # exact sampler with the time-average restriction mean(mu1*) < 0 (see eokim_exact_gibbs)
    ekc = do.call(eokim_gibbs, c(list(y, ndraw, burn, prior = ekp, seed = sd, verbose = FALSE, const_var = TRUE), EKA)))     # Eo-Kim paper-style sampler, constant error variance
  cols <- list()                                          # scalar parameter draws (keep rows, <= 6 columns)
  for (nm in setdiff(names(m), "state")) { v <- m[[nm]]
    if (is.numeric(v) && !is.null(dim(v)) && nrow(v) == keep && ncol(v) <= 6) { for (j in 1:ncol(v)) cols[[sprintf("%s%d", nm, j)]] <- v[, j] }
    else if (is.numeric(v) && is.null(dim(v)) && length(v) == keep) cols[[nm]] <- v }
  if (!is.null(cols$nu) && length(unique(cols$nu)) == 1) cols$nu <- NULL         # nu is fixed (5) in the t models: not a sampled parameter
  nmap <- c(alpha1 = "a0", alpha2 = "a1", P1 = "P00", P2 = "P11", theta1 = "theta1", theta2 = "theta2")
  D <- do.call(cbind, cols); cn <- colnames(D); cn[cn %in% names(nmap)] <- nmap[cn[cn %in% names(nmap)]]; colnames(D) <- cn
  if (!is.null(m$alpha)) D <- cbind(D, low = pmin(m$alpha[, 1], m$alpha[, 2]))[, , drop = FALSE]
  list(par = D, Sprob = m$Sprob, h = if (model$fam == "ekc") NULL else m$h, m0 = if (!is.null(m$m0)) mean(m$m0), m1 = if (!is.null(m$m1)) mean(m$m1), m0v = m$m0, m1v = m$m1, dv = m$delta,   # m0v/m1v: posterior-mean paths of the regime means (Eo-Kim)
      
       swap = m$swap_share, secs = proc.time()[3] - t0)
}

models <- list(); add <- function(nm, fam, tau = NA, rs = FALSE) models[[length(models) + 1]] <<- list(name = nm, fam = fam, tau = tau, restricted = rs)
add("Eo-Kim exact", "ekx"); add("Eo-Kim", "ek"); add("Eo-Kim const var", "ekc"); add("Eo-Kim exact (a1<0)", "ekxr", rs = TRUE)
for (rs in c(FALSE, TRUE)) { suf <- if (rs) " (a1<0)" else ""
  add(paste0("MSM-t-SV", suf), "msmtsv", rs = rs); add(paste0("MSM-SV", suf), "msmsv", rs = rs)
  add(paste0("MSM-t", suf), "msmt", rs = rs); add(paste0("MSM", suf), "msm", rs = rs)
  for (tau in seq(0.20, 0.80, by = 0.05)) add(sprintf("MSQ(%.2f)%s", tau, suf), "msq", tau, rs) }
if (drift) models <- Filter(function(m) m$fam %in% (if (tag == "driftq") "msq" else c("msm", "msmt", "msmsv", "msmtsv")) && !m$restricted, models)   # tag "driftq": the MSQ models with the slow drift (same zero-mean split)   # alpha1 < 0 is built into the zero-mean split
if (tag == "rwec") {      # the MSM-RW-EC family: base = constant variance, NO random-walk long-run mean (drift), NO sign restrictions; variants -t, -SV, -t-SV; toggles +drift (random-walk long-run mean) and (a1<0) (the paper's sign restrictions mu0* > 0 > mu1*)
  models <- list(); for (sv in c(FALSE, TRUE)) for (tt in c(FALSE, TRUE)) for (dr in c(FALSE, TRUE)) for (rs in c(FALSE, TRUE))
    models[[length(models) + 1]] <- list(name = paste0("MSM-RW-EC", if (tt) "-t" else "", if (sv) "-SV" else "", if (dr) " +drift" else "", if (rs) " (a1<0)" else ""), fam = "rwec", sv = sv, t = tt, drift = dr, restricted = rs, tau = NA) }
if (tag == "rep") {       # replication pages: Eo-Kim sample (vin = eokim) / Piger's factor (vin = vintage)
  models <- if (vin == "eokim") list(list(name = "Eo-Kim (paper sampler)", fam = "repek", restricted = FALSE, tau = NA), list(name = "Eo-Kim (paper sampler, our priors)", fam = "repek_ours", restricted = FALSE, tau = NA),
                                     list(name = "Eo-Kim exact (restricted)", fam = "repekx", restricted = FALSE, tau = NA), list(name = "MSM", fam = "msm", restricted = FALSE, tau = NA))
            else c(list(list(name = "MSM", fam = "msm", restricted = FALSE, tau = NA)), lapply(seq(0.20, 0.80, by = 0.05), function(t) list(name = sprintf("MSQ(%.2f)", t), fam = "msq", restricted = FALSE, tau = t))) }
if (nzchar(filt)) models <- Filter(function(m) grepl(filt, m$name), models)
tasks <- expand.grid(mi = seq_along(models), k = 1:K); cat(sprintf("[%s/%s] %d models x %d chains = %d tasks, %d draws (%d burn)\n", vin, tag, length(models), K, nrow(tasks), ndraw, burn))
t0 <- proc.time()[3]
res <- mclapply(seq_len(nrow(tasks)), function(i) { o <- tryCatch(fit1(models[[tasks$mi[i]]], tasks$k[i]), error = function(e) e)
  message(sprintf("[%4.0fs] %-18s chain %d %s", proc.time()[3] - t0, models[[tasks$mi[i]]]$name, tasks$k[i], if (inherits(o, "error")) paste("ERROR", conditionMessage(o)) else sprintf("(%.0fs)", o$secs))); o },
  mc.cores = 8, mc.preschedule = FALSE)
# post-processing (diagnostics, pooled results) is done per model, in parallel
post1 <- function(mi) { M <- models[[mi]]; fits <- res[tasks$mi == mi]; if (any(sapply(fits, inherits, "error"))) { return(list(error = TRUE)) }
  draws <- lapply(fits, `[[`, "par"); sp <- sapply(fits, `[[`, "Sprob"); pooled <- rowMeans(sp)
  W <- pmax(rowMeans(sapply(1:K, function(k) sp[, k] * (1 - sp[, k]))), 1e-12); Bv <- apply(sp, 1, var); rh <- sqrt(((keep - 1) / keep * W + Bv / keep) / W)
  summ <- chain_summary(draws); ps <- tryCatch(psrf_multi(lapply(draws, function(D) D[, setdiff(colnames(D), "low"), drop = FALSE])), error = function(e) NA)
  hh <- if (is.null(fits[[1]]$h)) rep(0, n) else rowMeans(sapply(fits, `[[`, "h")); a_pool <- if (!is.null(fits[[1]]$m0)) c(mean(sapply(fits, `[[`, "m0")), mean(sapply(fits, `[[`, "m1"))) else colMeans(do.call(rbind, lapply(draws, function(D) D[, c("a0", "a1")])))
  path <- if (!is.null(fits[[1]]$m0v)) list(m0 = rowMeans(sapply(fits, `[[`, "m0v")), m1 = rowMeans(sapply(fits, `[[`, "m1v"))) else NULL
  list(path = path, delta = if (!is.null(fits[[1]]$dv)) rowMeans(sapply(fits, `[[`, "dv")) else NULL, prob = pooled, alpha = as.numeric(a_pool), metrics = classify_metrics(pooled, S0), h = if (any(hh != 0)) hh else NULL,
    summary = summ, psrf = ps, S_rhat_max = max(rh), S_rhat_mean = mean(rh), S_maxdiff = max(apply(sp, 1, function(v) diff(range(v)))),
    chain_meanPr = colMeans(sp), chain_metrics = t(apply(sp, 2, classify_metrics, nber = S0)), chain_alpha = if (!is.null(fits[[1]]$m0)) NULL else t(sapply(draws, function(D) colMeans(D[, c("a0", "a1")]))),
    chain_Sprob = sp, secs = sapply(fits, `[[`, "secs"), swap = sapply(fits, function(f) if (is.null(f$swap)) NA else f$swap),
    nu = if ("nu" %in% colnames(draws[[1]])) mean(sapply(draws, function(D) mean(D[, "nu"]))))
}
out <- mclapply(seq_along(models), post1, mc.cores = 8, mc.preschedule = FALSE); names(out) <- sapply(models, `[[`, "name")
saveRDS(list(vintage = vin, tag = tag, K = K, ndraw = ndraw, burn = burn, date = r$date, S = S0, factor = y, models = out), sprintf("results/dfm_chains_%s_%s.rds", vin, tag))
cat(sprintf("\n%-20s %6s %6s %7s %7s | %5s %5s | %s\n", "model", "maxRhat", "PSRF", "minESSb", "minESSt", "S_Rhat", "dPr", "chain meanPr / AUC"))
for (nm in names(out)) { o <- out[[nm]]; if (isTRUE(o$error)) { cat(nm, "ERROR\n"); next }; s <- o$summary
  cat(sprintf("%-20s %6.2f %6.2f %7.0f %7.0f | %5.2f %5.2f | %s / %s\n", nm, max(s$rhat), o$psrf, min(s$ess_bulk, na.rm = TRUE), min(s$ess_tail, na.rm = TRUE), o$S_rhat_max, o$S_maxdiff,
      paste(round(o$chain_meanPr, 2), collapse = " "), paste(round(o$chain_metrics[, "AUC"], 2), collapse = " "))) }
cat(sprintf("total %.0fs\n", proc.time()[3] - t0))
