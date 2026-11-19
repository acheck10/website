# Univariate models on U.S. real GDP growth (FRED GDPC1), scored against NBER recession quarters (FRED USRECQ).
#
# Growth = compounded annualized rate of change: 100 * ((GDP_t / GDP_{t-1})^4 - 1)   (FRED "pca").
# Scaling: all priors (regime means N(+/-0.5, 1), Eo-Kim variance priors, ...) are set for QUARTERLY percent growth,
# so the models are fit to growth / 4 and mean-type parameters are multiplied by 4 on the way out. Every model here
# is scale-equivariant (improper sigma priors, data-based SV starting values), so this is identical to rescaling
# every prior by 4.
#
# usage: Rscript scripts/run_gdp.R [ndraw] [burn]     -> gdp_results.rds, gdp_results.json
args <- commandArgs(trailingOnly = TRUE)
ndraw <- if (length(args) >= 1) as.integer(args[1]) else 6000L
burn  <- if (length(args) >= 2) as.integer(args[2]) else 2000L
if (nzchar(Sys.getenv("PULL_FRED"))) source("scripts/pull_fred.R")   # PULL_FRED=1: refresh the "current" data (latest vintage) with fredr first
source("R/ms_samplers.R"); source("R/ms_models_extra.R"); source("R/ms_models_sv.R"); source("R/eokim_exact.R"); source("R/ms_models_sv2.R"); source("R/ms_models_drift0.R")
SC <- 4                                                     # annualisation factor
# A1MAX=0 (environment variable) runs the univariate models with the restriction a1 < 0 on the low-regime mean
# (in fitted units, i.e. growth/4; 0 is 0 either way) and skips the Eo-Kim samplers (their own sign restrictions apply).
A1MAX <- Sys.getenv("A1MAX", ""); restricted <- nzchar(A1MAX)
# DRIFT=1 (environment variable) fits only the four switching-mean models with the slow drift delta_t (Eo-Kim's random-walk long-run mean, fixed variance 0.0004)
# and the zero-mean split pi0*a0 + pi1*a1 = 0 (R/ms_models_drift0.R; a1 ~ N(-1.5, 1) truncated < 0). A1MAX is not used with DRIFT. Outputs gdp_results_drift.{rds,json}.
drift <- nzchar(Sys.getenv("DRIFT", ""))
# CONSTVAR=1 fits only Eo-Kim's paper-style sampler with a CONSTANT error variance instead of stochastic volatility -> gdp_results_constvar.{rds,json}
constvar <- nzchar(Sys.getenv("CONSTVAR", ""))
pr  <- modifyList(default_prior, if (restricted) list(a1_max = as.numeric(A1MAX)) else list())
prs <- modifyList(sv_prior,      if (restricted) list(a1_max = as.numeric(A1MAX)) else list())

gdp <- read.csv("data_gdpc1/GDPC1.csv"); rec <- read.csv("data_usrecq/USRECQ.csv")
names(gdp) <- c("date", "gdp"); names(rec) <- c("date", "rec")
d <- merge(gdp, rec, by = "date")
d$g <- c(NA, 100 * ((d$gdp[-1] / d$gdp[-nrow(d)])^4 - 1))   # compounded annualised growth
d <- d[!is.na(d$g), ]
d$date <- as.Date(d$date)
cat(sprintf("Sample %s to %s, T = %d quarters; %d NBER recession quarters. Growth: mean %.2f, sd %.2f, min %.1f (%s)\n",
            format(d$date[1], "%Y-Q%q"), format(tail(d$date, 1)), nrow(d), sum(d$rec), mean(d$g), sd(d$g),
            min(d$g), format(d$date[which.min(d$g)])))

fit_all <- function(dd, label) {
  y <- dd$g / SC; nber <- dd$rec; out <- list(); t0 <- proc.time()[3]
  tic <- function(nm) message(sprintf("[%s] %-10s done (%.0fs)", label, nm, proc.time()[3] - t0))
  rep_m <- function(prob, alpha4, extra = list()) c(list(prob = prob, alpha = alpha4,
                    metrics = classify_metrics(prob, nber)), extra)
  if (constvar) {
    m <- eokim_gibbs(y, ndraw, burn, const_var = TRUE, seed = 1, verbose = FALSE); tic("Eo-Kim const var")
    out$`Eo-Kim const var` <- rep_m(m$Sprob, SC * c(mean(m$m0), mean(m$m1)), list(path = list(m0 = SC * m$m0, m1 = SC * m$m1), theta = colMeans(m$theta), s2w = colMeans(m$s2w), P = colMeans(m$P)))
    return(list(date = format(dd$date), growth = dd$g, nber = nber, models = out)) }
  if (drift) {
    for (nm in c("MSM", "MSM-t", "MSM-SV", "MSM-t-SV")) { fam <- c(MSM = "msm", `MSM-t` = "msmt", `MSM-SV` = "msmsv", `MSM-t-SV` = "msmtsv")[[nm]]
      m <- drift0_gibbs(fam, y, ndraw, burn, prior = if (fam %in% c("msm", "msmt")) pr else prs, a1_prior = c(-1.5, 1), s2d = 4e-4, d0 = mean(y), dsd = 0.5, seed = 1, verbose = FALSE)
      out[[nm]] <- rep_m(m$Sprob, SC * colMeans(m$alpha), c(list(delta = SC * m$delta, P = colMeans(m$P)), if (!is.null(m$nu)) list(nu = mean(m$nu)), if (!is.null(m$h)) list(sd = SC * exp(m$h / 2))))
      tic(paste(nm, "+drift")) }
    return(list(date = format(dd$date), growth = dd$g, nber = nber, models = out)) }
  m <- msm_gibbs(y, ndraw, burn, prior = pr, seed = 1, verbose = FALSE);   tic("MSM")
  out$MSM <- rep_m(m$Sprob, SC * colMeans(m$alpha), list(sigma = SC * mean(m$sigma), P = colMeans(m$P)))
  m <- msmt_gibbs(y, ndraw, burn, prior = pr, seed = 1, verbose = FALSE);  tic("MSM-t")
  out$`MSM-t` <- rep_m(m$Sprob, SC * colMeans(m$alpha), list(nu = mean(m$nu), P = colMeans(m$P)))
  m <- msmsv_gibbs(y, ndraw, burn, prior = prs, seed = 1, verbose = FALSE); tic("MSM-SV")
  out$`MSM-SV` <- rep_m(m$Sprob, SC * colMeans(m$alpha), list(sd = SC * exp(m$h / 2), P = colMeans(m$P)))
  m <- msmtsv_gibbs(y, ndraw, burn, prior = prs, seed = 1, verbose = FALSE); tic("MSM-t-SV")
  out$`MSM-t-SV` <- rep_m(m$Sprob, SC * colMeans(m$alpha), list(nu = mean(m$nu), sd = SC * exp(m$h / 2), P = colMeans(m$P)))
  for (tau in seq(0.20, 0.80, by = 0.05)) {
    m <- msq_gibbs(y, tau = tau, ndraw, burn, prior = pr, seed = 1, verbose = FALSE)
    out[[sprintf("MSQ(%.2f)", tau)]] <- rep_m(m$Sprob, SC * colMeans(m$alpha), list(P = colMeans(m$P)))
  }; tic("MSQ grid")
  if (restricted) return(list(date = format(dd$date), growth = dd$g, nber = nber, models = out))
  m <- eokim_gibbs(y, ndraw, burn, seed = 1, verbose = FALSE); tic("Eo-Kim")
  out$`Eo-Kim` <- rep_m(m$Sprob, SC * c(mean(m$m0), mean(m$m1)),
                        list(path = list(m0 = SC * m$m0, m1 = SC * m$m1), theta = colMeans(m$theta),
                             s2w = colMeans(m$s2w), P = colMeans(m$P)))
  m <- eokim_exact_gibbs(y, ndraw, burn, seed = 1, verbose = FALSE); tic("Eo-Kim exact")
  out$`Eo-Kim exact` <- rep_m(m$Sprob, SC * c(mean(m$m0), mean(m$m1)),
                        list(path = list(m0 = SC * m$m0, m1 = SC * m$m1), theta = colMeans(m$theta),
                             s2w = colMeans(m$s2w), P = colMeans(m$P), swap_share = m$swap_share, acc = m$acc))
  list(date = format(dd$date), growth = dd$g, nber = nber, models = out)
}

res <- list(full = fit_all(d, "1947Q2-"))
# the Eo-Kim paper's sample, 1947Q4-2011Q3, as a replication check
d2 <- d[d$date >= as.Date("1947-10-01") & d$date <= as.Date("2011-07-01"), ]
res$eokim_sample <- fit_all(d2, "1947Q4-2011Q3")
sfx <- paste0(if (drift) "_drift" else "", if (constvar) "_constvar" else "", if (restricted) "_a1neg" else "")
saveRDS(res, sprintf("results/gdp_results%s.rds", sfx))
writeLines(jsonlite::toJSON(res, auto_unbox = TRUE, digits = 6), sprintf("webpage/gdp_results%s.json", sfx))

tab <- function(r) do.call(rbind, lapply(names(r$models), function(n) { m <- r$models[[n]]$metrics
  data.frame(model = n, AUC = m["AUC"], TPR = m["TPR"], FPR = m["FPR"], Prec = m["Precision"], F1 = m["F1"],
             low = r$models[[n]]$alpha[2], high = r$models[[n]]$alpha[1], row.names = NULL) }))
for (nm in names(res)) { cat("\n====", nm, "====\n"); print(format(tab(res[[nm]]), digits = 2), row.names = FALSE) }
