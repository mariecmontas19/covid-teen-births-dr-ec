# 09_monthly_conceptions.R — short-term Kearney & Levine (2023)-style monthly
# conception analysis for the DR, mirroring Omar's Ecuador draft: fit
# log(monthly conceptions) = month-of-year FE + linear trend ON PRE-PANDEMIC
# DATA ONLY (2016m1-2020m2), then compute out-of-sample deviations of actual
# from predicted conceptions for 2020m3-2024m12. Conceptions dated via the
# microdata conception_date. Window ends 2024m12 (fully observed).
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

b <- as.data.table(readRDS(file.path(DR_CLEAN, "births_clean_2016_2025.rds")))
b <- b[!is.na(conception_date) & !is.na(adm3_pcode)]
b[, `:=`(cy = year(conception_date), cm = month(conception_date))]
b <- b[cy >= 2016 & cy <= 2024]
b[, teen := age_mom >= 15 & age_mom <= 19]

w10 <- as.data.table(readRDS(file.path(OUT, "wealth2010_muni.rds")))
q20 <- quantile(w10$wealth_index_2010, .20)
w10[, poorest := as.integer(wealth_index_2010 <= q20)]
b <- merge(b, w10[, .(adm3_pcode, poorest)], by = "adm3_pcode", all.x = TRUE)

# K&L procedure: pre-period model, out-of-sample deviations.
# DATA DECISION (MM, 2026-07-27): the registry ramps 104k (2016) -> 162k (2019)
# births because the MOH was phasing in the new live-birth registration system.
# The monthly analysis therefore uses 2019m1-2020m2 as the pre window (post-
# phase-in), with month fixed effects ONLY (14 months cannot support a trend),
# so deviations are relative to 2019 seasonal means. Horizon: through 2021m12.
# The poorest-vs-rest DIFFERENCE in deviations is the credible object here —
# registration frictions common to both groups cancel.
PRE_START <- as.Date("2019-01-01")
kl_dev <- function(dt) {
  d <- dt[, .(n = .N), by = .(cy, cm)][order(cy, cm)]
  d[, `:=`(t = (cy - 2016) * 12 + cm, logn = log(n),
           date = as.Date(sprintf("%d-%02d-01", cy, cm)))]
  pre <- d[date >= PRE_START & date < as.Date("2020-03-01")]
  m <- lm(logn ~ factor(cm), data = pre)   # month FE only; no trend
  d[, pred := predict(m, newdata = d)]
  d[, dev := logn - pred]
  d[, se_pre := summary(m)$sigma]           # pre-period residual SD as band
  d
}

samples <- list(teen_all   = b[teen == TRUE],
                teen_poor  = b[teen == TRUE & poorest == 1],
                teen_rich  = b[teen == TRUE & poorest == 0],
                adult_2029 = b[age_mom >= 20 & age_mom <= 29])
devs <- rbindlist(lapply(names(samples), function(nm) kl_dev(samples[[nm]])[, grp := nm]))
saveRDS(devs, file.path(OUT, "dr_monthly_kl.rds"))

# Summary: missing/excess conceptions in key windows
sink(file.path(OUT, "logs", "09_monthly_conceptions.txt"))
for (nm in names(samples)) {
  d <- devs[grp == nm]
  win <- function(a, z) {
    x <- d[date >= a & date <= z]
    sprintf("%+.1f%% (%+d births)", 100 * (exp(mean(x$dev)) - 1), round(sum(exp(x$pred) * (exp(x$dev) - 1))))
  }
  cat(sprintf("%-11s bust 2020m3-2020m12: %s | 2021: %s\n", nm,
              win("2020-03-01", "2020-12-01"), win("2021-01-01", "2021-12-01")))
}
# Poorest-vs-rest differential in deviations (registration frictions cancel)
dd <- dcast(devs[grp %in% c("teen_poor", "teen_rich"), .(grp, date, dev)],
            date ~ grp, value.var = "dev")[, dif := teen_poor - teen_rich]
cat(sprintf("\nPoor-rich teen differential: 2020m3-m12 %+.1f pp | 2021 %+.1f pp\n",
            100 * dd[date >= as.Date("2020-03-01") & date <= as.Date("2020-12-01"), mean(dif)],
            100 * dd[date >= as.Date("2021-01-01") & date <= as.Date("2021-12-01"), mean(dif)]))
sink()
system(paste("cat", shQuote(file.path(OUT, "logs", "09_monthly_conceptions.txt"))))

# Table 5 (replaces the earlier monthly figure — MM 2026-07-27: figure too busy,
# present the windows as a compact table instead): deviation of conceptions from
# the 2019 seasonal baseline, by group and window, with implied birth counts.
cell <- function(nm, a, z) {
  x <- devs[grp == nm & date >= a & date <= z]
  c(pct = 100 * (exp(mean(x$dev)) - 1), n = sum(exp(x$pred) * (exp(x$dev) - 1)))
}
rows <- list(
  c("Teens 15--19 (all)",          "teen_all"),
  c("\\quad Poorest 20\\%",        "teen_poor"),
  c("\\quad Other 80\\%",          "teen_rich"),
  c("Women 20--29",                "adult_2029"))
lines <- vapply(rows, function(r) {
  b <- cell(r[2], as.Date("2020-03-01"), as.Date("2020-12-01"))
  y <- cell(r[2], as.Date("2021-01-01"), as.Date("2021-12-01"))
  sprintf("%s & %+.1f\\%% & %s & %+.1f\\%% & %s \\\\", r[1],
          b["pct"], format(round(b["n"]), big.mark = ","),
          y["pct"], format(round(y["n"]), big.mark = ","))
}, character(1))
tab <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Dominican Republic: conceptions relative to the 2019 seasonal baseline}",
  "\\label{tab:dr_monthly}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Mar--Dec 2020 conceptions} & \\multicolumn{2}{c}{2021 conceptions} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  " & \\% deviation & births & \\% deviation & births \\\\",
  "\\midrule",
  lines,
  "\\midrule",
  sprintf("Poorest $-$ other teen differential & \\multicolumn{2}{c}{%+.1f pp} & \\multicolumn{2}{c}{%+.1f pp} \\\\",
          100 * dd[date >= as.Date("2020-03-01") & date <= as.Date("2020-12-01"), mean(dif)],
          100 * dd[date >= as.Date("2021-01-01") & date <= as.Date("2021-12-01"), mean(dif)]),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.92\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Conceptions dated from birth records (birth date minus gestation).",
  "Baseline: month-of-year means estimated on 2019m1--2020m2, after the",
  "registry phase-in of 2016--2019. Deviations are out of sample; birth counts",
  "are implied excess/missing conceptions. Because no trend is removed, level",
  "deviations partly reflect the secular decline in teen fertility; the",
  "poorest-vs-other differential nets out registry and secular components",
  "common to both groups.",
  "\\end{minipage}",
  "\\end{table}")
writeLines(tab, file.path(OUT_TAB, "tab5_dr_monthly_conceptions.tex"))
cat("Saved tab5_dr_monthly_conceptions.tex and summary log.\n")
