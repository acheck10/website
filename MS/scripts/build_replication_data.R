# Eo & Kim replication dataset: real GDP growth exactly as the authors' own data file defines it, from the vintage they would have used.
#   growth = simple quarterly percent change of GDPC1, ALFRED vintage 2011-12-22 (the third estimate of 2011Q3; data_gdpc1_eokim/GDPC1_alfred_2011-12-22.csv, pulled with fredr): 1947Q4 - 2011Q3, 256 quarters
#   S      = the authors' recession indicator (eo_kim_replication_code/realgdp11q3_recess.txt: includes the peak quarter, drops the trough quarter; differs from FRED's USRECQ in 8 quarters)
# The ALFRED series reproduces the authors' growth file to within 0.004 percentage points (median 0.0004); their file is kept as growth_authors for reference.
# Output: results/eokim_data.rds in the format run_dfm_chains.R expects (factor = growth, fitted as is: the series is already on the scale the paper's priors were written for).   usage: Rscript scripts/build_replication_data.R   (OUTDIR=<dir>; DRY=1 prints only)
out <- Sys.getenv("OUTDIR", "results")
v <- read.csv("data_gdpc1_eokim/GDPC1_alfred_2011-12-22.csv"); v$observation_date <- as.Date(v$observation_date); v <- v[!is.na(v$GDPC1), ]
g <- 100 * (v$GDPC1[-1] / v$GDPC1[-nrow(v)] - 1); d <- v$observation_date[-1]; i <- d >= as.Date("1947-10-01") & d <= as.Date("2011-07-01"); g <- g[i]; d <- d[i]
a <- read.table("eo_kim_replication_code/realgdp11q3_recess.txt"); names(a) <- c("g", "rec"); stopifnot(length(g) == nrow(a))
cat(sprintf("eokim: %d quarters %s to %s; growth mean %.3f sd %.3f; recession quarters (authors) %d; max |ALFRED - authors| = %.4f\n", length(g), format(d[1]), format(tail(d, 1)), mean(g), sd(g), sum(a$rec), max(abs(g - a$g))))
if (!nzchar(Sys.getenv("DRY"))) saveRDS(list(vintage = "eokim", date = format(d), S = as.integer(a$rec), growth = g, factor = g, growth_authors = a$g, alfred_vintage = "2011-12-22"), file.path(out, "eokim_data.rds"))
