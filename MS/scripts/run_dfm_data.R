# Two-step dynamic-factor models on the four series in Piger (2026): PAYEMS, INDPRO ("IP"), W875RX1, CMRMTSPL, scored against NBER (USREC).
# usage: Rscript scripts/run_dfm_data.R <vintage|current> [ndraw] [burn]     vintage = ALFRED 2025-09-30 (the paper's final vintage), Feb 1967 - Jul 2025
a <- commandArgs(trailingOnly = TRUE); vin <- a[1]; nd <- if (length(a) >= 2) as.integer(a[2]) else 3000L; bn <- if (length(a) >= 3) as.integer(a[3]) else 1000L
if (vin == "current" && nzchar(Sys.getenv("PULL_FRED"))) source("scripts/pull_fred.R")   # PULL_FRED=1: refresh the current data (latest vintage) with fredr first
source("R/ms_samplers.R"); source("R/ms_models_extra.R"); source("R/ms_models_sv.R"); source("R/eokim_exact.R"); source("R/df_factor.R")
ids <- c("PAYEMS", "INDPRO", "W875RX1", "CMRMTSPL")
rd <- function(id) { f <- if (vin == "vintage") sprintf("data_dfm/%s_vintage/%s_2025-09-30.csv", id, id) else sprintf("data_dfm/%s/%s.csv", id, id)
  d <- read.csv(f, stringsAsFactors = FALSE); names(d) <- c("date", "v"); d$v <- suppressWarnings(as.numeric(d$v)); d[!is.na(d$v), ] }
gr <- function(d) data.frame(date = d$date[-1], g = 100 * diff(log(d$v)))                    # monthly growth, 100 x log difference
G <- Reduce(function(x, y) merge(x, y, by = "date"), lapply(ids, function(i) { g <- gr(rd(i)); names(g)[2] <- i; g }))
usrec <- read.csv("data_dfm/USREC/USREC.csv"); names(usrec) <- c("date", "rec"); G <- merge(G, usrec, by = "date")
if (vin == "vintage") G <- G[G$date <= "2025-07-01", ]
n <- nrow(G); S <- G$rec; Xraw <- as.matrix(G[, ids]); X <- scale(Xraw)                      # demeaned, unit-variance series for the factor step
cat(sprintf("[%s] %d months, %s to %s; NBER recession months %d (%.0f%%); series means (monthly %%): %s\n", vin, n, G$date[1], G$date[n], sum(S), 100 * mean(S), paste(round(colMeans(Xraw), 2), collapse = " ")))
t1 <- proc.time()[3]; tic <- function(m) message(sprintf("[%s] %-14s done (%.0fs)", vin, m, proc.time()[3] - t1))
fac <- dfactor_gibbs(X, nd, bn, seed = 1); y <- fac$Ymed; tic("factor")
ph <- colMeans(fac$phi); cat(sprintf("factor step: loadings %s | idio var %s | phi %s (phi1+phi2 = %.2f; bound %.2f) | corr(factor, NBER indicator) = %.2f\n",
    paste(round(colMeans(fac$psi), 2), collapse = ","), paste(round(colMeans(fac$sig2), 2), collapse = ","), paste(round(ph, 2), collapse = ","), sum(ph), df_prior$rho, cor(y, S)))
cat(sprintf("factor (posterior median): sd %.2f | min %.1f at %s | mean in NBER recessions %.2f, in expansions %.2f\n", sd(y), min(y), G$date[which.min(y)], mean(y[S == 1]), mean(y[S == 0])))
res <- list(); add <- function(nm, prob, alpha, extra = list()) res[[nm]] <<- c(list(prob = prob, alpha = alpha, metrics = classify_metrics(prob, S)), extra)
rest <- list(a1_max = 0)
for (rs in c(FALSE, TRUE)) { suf <- if (rs) " (a1<0)" else ""; pr <- modifyList(default_prior, if (rs) rest else list()); prs <- modifyList(sv_prior, if (rs) rest else list())
  m <- msm_gibbs(y, nd, bn, prior = pr, seed = 1, verbose = FALSE);   add(paste0("MSM", suf), m$Sprob, colMeans(m$alpha))
  m <- msmt_gibbs(y, nd, bn, prior = pr, seed = 1, verbose = FALSE);  add(paste0("MSM-t", suf), m$Sprob, colMeans(m$alpha), list(nu = mean(m$nu)))
  m <- msmsv_gibbs(y, nd, bn, prior = prs, seed = 1, verbose = FALSE); add(paste0("MSM-SV", suf), m$Sprob, colMeans(m$alpha))
  m <- msmtsv_gibbs(y, nd, bn, prior = prs, seed = 1, verbose = FALSE); add(paste0("MSM-t-SV", suf), m$Sprob, colMeans(m$alpha), list(nu = mean(m$nu)))
  for (tau in seq(0.20, 0.80, by = 0.05)) { m <- msq_gibbs(y, tau = tau, nd, bn, prior = pr, seed = 1, verbose = FALSE); add(sprintf("MSQ(%.2f)%s", tau, suf), m$Sprob, colMeans(m$alpha)) }
  tic(if (rs) "restricted" else "unrestricted") }
m <- eokim_gibbs(y, nd, bn, seed = 1, verbose = FALSE); add("Eo-Kim", m$Sprob, c(mean(m$m0), mean(m$m1)), list(theta = colMeans(m$theta))); tic("Eo-Kim")
m <- eokim_exact_gibbs(y, nd, bn, seed = 1, verbose = FALSE); add("Eo-Kim exact", m$Sprob, c(mean(m$m0), mean(m$m1)), list(swap_share = m$swap_share)); tic("Eo-Kim exact")
saveRDS(list(vintage = vin, date = G$date, S = S, X = Xraw, factor = y, models = res), sprintf("results/dfm_data_%s.rds", vin))
tab <- do.call(rbind, lapply(names(res), function(k) { q <- res[[k]]; data.frame(model = k, AUC = q$metrics["AUC"], TPR = q$metrics["TPR"], FPR = q$metrics["FPR"], Prec = q$metrics["Precision"], F1 = q$metrics["F1"], low = q$alpha[2], high = q$alpha[1], row.names = NULL) }))
cat("\n"); print(format(tab, digits = 2), row.names = FALSE)
