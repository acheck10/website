# Builds the stress-test datasets for the four-series factor page (all start from results/dfm_data_current.rds: monthly growth of PAYEMS, INDPRO, W875RX1, CMRMTSPL, Feb 1967 - Aug 2026).
#   cens : standard one-factor model; the five COVID months (Mar-Jul 2020) are CENSORED (dropped) before the regime models are fitted.
#   adj  : factor model with a parameter change for Mar-Jul 2020 (R/df_factor_covid.R, after Piger's FAQ) then the COVID months are dropped for the regime models.
#   e90c / e01c : the same extensions with the standard factor and the COVID months censored (as "cens").
#   e90 / e01 : the real data EXTENDED with a simulated mild recession starting in the first month out of sample (Oct 2026), shaped like 1990-91 / 2001, standard factor model, COVID months kept
#               (this is what the pages' models face today, plus the new recession);  e90a / e01a: the same extensions with the "adj" treatment.
# Simulated extension: for each series, ext_t = recent mean (last 36 months) + (template_t - template pre-recession mean) + N(0, (0.5 s_i)^2), where the template is the actual 1990 / 2001 episode
# (one expansion month, the recession months, the recovery and then ONE MORE YEAR OF EXPANSION: 42 months in all, Sep 2026 - Feb 2030; EXT_H overrides, 30 = the earlier version without the extra year) and s_i is the series' month-to-month sd in post-1984 expansions: a real recession pattern on today's trend growth, plus extra noise.
# usage: Rscript scripts/build_stress_data.R        -> results/stress_<name>_data.rds
source("R/df_factor.R"); source("R/df_factor_covid.R")
r <- readRDS("results/dfm_data_current.rds"); X0 <- r$X; d0 <- as.Date(r$date); S0 <- r$S; n0 <- nrow(X0); N <- ncol(X0); colnames(X0) <- c("PAYEMS", "INDPRO", "W875RX1", "CMRMTSPL")
covid_dates <- seq(as.Date("2020-03-01"), as.Date("2020-07-01"), by = "month")
H_EXT <- as.integer(Sys.getenv("EXT_H", "42"))
ext_X <- function(tpl_first_rec, seed, H = H_EXT) {
  set.seed(seed); r0 <- which(d0 == as.Date(tpl_first_rec)); stopifnot(S0[r0] == 1, S0[r0 - 1] == 0)
  rows <- (r0 - 1):(r0 - 2 + H); tpl <- X0[rows, , drop = FALSE]; pre <- X0[(r0 - 61):(r0 - 2), , drop = FALSE]
  recent <- colMeans(X0[(n0 - 35):n0, ]); calm <- d0 >= as.Date("1985-01-01") & S0 == 0 & !(d0 %in% covid_dates); s <- apply(X0[calm, ], 2, sd)
  Xe <- t(t(tpl) - colMeans(pre)) ; Xe <- t(t(Xe) + recent) + matrix(rnorm(H * N, 0, rep(0.5 * s, each = H)), H, N)
  list(X = Xe, S = S0[rows], date = seq(d0[n0] + 32, by = "month", length.out = H) - as.integer(format(seq(d0[n0] + 32, by = "month", length.out = H), "%d")) + 1) }
ext <- list(e90 = ext_X("1990-08-01", 901), e01 = ext_X("2001-04-01", 2001))
if (nzchar(Sys.getenv("DRY"))) { for (nm in names(ext)) { e <- ext[[nm]]; cat(sprintf("[DRY] %s: %d months %s to %s; NBER months %d (first %s); last %d months all expansion: %s\n", nm, length(e$S), format(e$date[1]), format(tail(e$date, 1)), sum(e$S), format(e$date[which(e$S == 1)[1]]), 12, all(tail(e$S, 12) == 0))) }; quit(save = "no") }
make <- function(name, X, S, d, adj, drop_covid, ext_from = NA) {
  Z <- scale(X); cv <- which(d %in% covid_dates)
  if (adj) { f <- dfactor_gibbs_cv(Z, 3000, 1000, covid = cv, seed = 1); y <- f$Ymed; o <- setdiff(seq_along(y), cv); y <- y * (0.63 / sd(y[o])); cat(sprintf("%-5s kappa %.0f, vc %.0f, phi %s\n", name, mean(f$kappa), mean(f$vc), paste(round(colMeans(f$phi), 2), collapse = ","))) }
  else { f <- dfactor_gibbs(Z, 3000, 1000, seed = 1); y <- f$Ymed; cat(sprintf("%-5s phi %s\n", name, paste(round(colMeans(f$phi), 2), collapse = ","))) }
  keep <- if (drop_covid) setdiff(seq_along(y), cv) else seq_along(y)
  saveRDS(list(vintage = paste0("stress_", name), date = format(d[keep]), S = S[keep], growth = y[keep], factor = y[keep], factor_full = y, dropped = if (drop_covid) format(d[cv]) else NULL, ext_from = ext_from), sprintf("results/stress_%s_data.rds", name))
  cat(sprintf("      T = %d, NBER months %d, factor sd %.2f, min %.1f\n", length(keep), sum(S[keep]), sd(y[keep]), min(y[keep]))) }
skip_adj <- nzchar(Sys.getenv("SKIP_ADJ"))      # SKIP_ADJ=1: do not build the COVID-adjusted samples (no longer shown on the page; each needs a factor fit with a parameter change)
make("cens", X0, S0, d0, FALSE, TRUE); if (!skip_adj) make("adj", X0, S0, d0, TRUE, TRUE)
for (nm in names(ext)) { e <- ext[[nm]]; Xa <- rbind(X0, e$X); Sa <- c(S0, e$S); da <- c(d0, e$date); fr <- format(e$date[1])
  cat(sprintf("%s extension %s to %s: %d NBER months from %s\n", nm, format(e$date[1]), format(tail(e$date, 1)), sum(e$S), format(e$date[which(e$S == 1)[1]])))
  make(nm, Xa, Sa, da, FALSE, FALSE, fr); if (!skip_adj) make(paste0(nm, "a"), Xa, Sa, da, TRUE, TRUE, fr); make(paste0(nm, "c"), Xa, Sa, da, FALSE, TRUE, fr) }
saveRDS(ext, "results/stress_extensions.rds")
