# 10_monthly_conceptions_both.R — harmonized short-run Kearney & Levine (2023)
# conception analysis for BOTH countries (approved MM/OG 2026-07-28).
# Conception dating: birth month minus 9 in both countries (OG's convention);
# the DR gestational-age dating remains in 09_monthly_conceptions.R as
# robustness. Pre-period models:
#   Ecuador:  month FE + linear trend, fit on 2010m1-2020m2 (spec A, OG draft)
#             and on 2016m1-2020m2 (spec B, harmonized with the DR frame).
#   DR:       month FE only, fit on 2019m1-2020m2 (registry phase-in 2016-19).
# Groups: teens 15-19 (all / poorest quintile / other), women 20-29; Ecuador
# adds all-ages to validate against OG's national bust (-13%, 38,472 missing).
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

shift9 <- function(dt) {           # conception ym = birth ym - 9 months
  dt[, cm := bm - 9L]
  dt[, cy := fifelse(cm <= 0L, by_ - 1L, by_)]
  dt[, cm := fifelse(cm <= 0L, cm + 12L, cm)]
  dt
}
kl_dev <- function(dt, pre_start, pre_end = as.Date("2020-02-01"), trend = TRUE) {
  d <- dt[, .(n = .N), by = .(cy, cm)][order(cy, cm)]
  d <- d[!is.na(cy)]
  d[, `:=`(t = (cy - min(cy)) * 12 + cm, logn = log(n),
           date = as.Date(sprintf("%d-%02d-01", cy, cm)))]
  pre <- d[date >= pre_start & date <= pre_end]
  f <- if (trend) logn ~ factor(cm) + t else logn ~ factor(cm)
  m <- lm(f, data = pre)
  d[, pred := predict(m, newdata = d)]
  d[, dev := logn - pred]
  d
}
win <- function(d, a, z) {
  x <- d[date >= a & date <= z]
  c(pct = 100 * (exp(mean(x$dev)) - 1), n = sum(exp(x$pred) * (exp(x$dev) - 1)))
}
W2020 <- c(as.Date("2020-03-01"), as.Date("2020-12-01"))
W2021 <- c(as.Date("2021-01-01"), as.Date("2021-12-01"))

# ---- Ecuador -----------------------------------------------------------------
ec <- as.data.table(haven::read_dta(file.path(EC_DATA, "ENV", "ENV_CLEAN.dta"),
                                    col_select = c("year", "month", "age", "canton")))
setnames(ec, c("year", "month"), c("by_", "bm"))
ec <- shift9(ec)[cy >= 2009 & cy <= 2023]
ec[, canton := as.numeric(canton)]
ai <- as.data.table(haven::read_dta(file.path(EC_DATA, "CPV2010", "CPV2010_aindex_canton.dta")))
ec <- merge(ec, ai[, .(canton = as.numeric(canton), poorest)], by = "canton", all.x = TRUE)

ec_samples <- list(
  all_ages  = ec,
  teen_all  = ec[age >= 15 & age <= 19],
  teen_poor = ec[age >= 15 & age <= 19 & poorest == 1],
  teen_rich = ec[age >= 15 & age <= 19 & poorest == 0],
  adult2029 = ec[age >= 20 & age <= 29])

ecA <- lapply(ec_samples, kl_dev, pre_start = as.Date("2010-01-01"), trend = TRUE)
ecB <- lapply(ec_samples, kl_dev, pre_start = as.Date("2016-01-01"), trend = TRUE)

# Validation vs OG draft: 2020 full-year missing conceptions, all ages, spec A
v <- ecA$all_ages[cy == 2020]
cat(sprintf("EC validation (spec A, all ages, calendar 2020): %+.1f%% (%+d conceptions) [OG: -13%%, -38,472]\n",
            100 * (exp(mean(v$dev)) - 1), round(sum(exp(v$pred) * (exp(v$dev) - 1)))))

# ---- Dominican Republic (9-month dating, harmonized) --------------------------
b <- as.data.table(readRDS(file.path(DR_CLEAN, "births_clean_2016_2025.rds")))
b <- b[!is.na(adm3_pcode) & !is.na(birth_year) & !is.na(birth_month)]
setnames(b, c("birth_year", "birth_month"), c("by_", "bm"))
b <- shift9(b)[cy >= 2016 & cy <= 2024]
w10 <- as.data.table(readRDS(file.path(OUT, "wealth2010_muni.rds")))
w10[, poorest := as.integer(wealth_index_2010 <= quantile(wealth_index_2010, .20))]
b <- merge(b, w10[, .(adm3_pcode, poorest)], by = "adm3_pcode", all.x = TRUE)

dr_samples <- list(
  all_ages  = b,
  teen_all  = b[age_mom >= 15 & age_mom <= 19],
  teen_poor = b[age_mom >= 15 & age_mom <= 19 & poorest == 1],
  teen_rich = b[age_mom >= 15 & age_mom <= 19 & poorest == 0],
  adult2029 = b[age_mom >= 20 & age_mom <= 29])
dr <- lapply(dr_samples, kl_dev, pre_start = as.Date("2019-01-01"), trend = FALSE)

# ---- Summary log ---------------------------------------------------------------
sink(file.path(OUT, "logs", "10_monthly_conceptions_both.txt"))
fmt <- function(w) sprintf("%+.1f%% (%+d)", w["pct"], round(w["n"]))
cat("== ECUADOR spec A (pre 2010m1-2020m2, month FE + trend) ==\n")
for (nm in names(ecA)) cat(sprintf("%-10s 2020m3-m12: %s | 2021: %s\n", nm,
    fmt(win(ecA[[nm]], W2020[1], W2020[2])), fmt(win(ecA[[nm]], W2021[1], W2021[2]))))
cat("\n== ECUADOR spec B (pre 2016m1-2020m2, month FE + trend) ==\n")
for (nm in names(ecB)) cat(sprintf("%-10s 2020m3-m12: %s | 2021: %s\n", nm,
    fmt(win(ecB[[nm]], W2020[1], W2020[2])), fmt(win(ecB[[nm]], W2021[1], W2021[2]))))
cat("\n== DR (9-month dating, pre 2019m1-2020m2, month FE only) ==\n")
for (nm in names(dr)) cat(sprintf("%-10s 2020m3-m12: %s | 2021: %s\n", nm,
    fmt(win(dr[[nm]], W2020[1], W2020[2])), fmt(win(dr[[nm]], W2021[1], W2021[2]))))
diff <- function(lst) c(
  100 * (win(lst$teen_poor, W2020[1], W2020[2])["pct"] - win(lst$teen_rich, W2020[1], W2020[2])["pct"]) / 100,
  100 * (win(lst$teen_poor, W2021[1], W2021[2])["pct"] - win(lst$teen_rich, W2021[1], W2021[2])["pct"]) / 100)
cat(sprintf("\nPoor-rich teen differential (pp): EC specB 2020 %+.1f | 2021 %+.1f ;  DR 2020 %+.1f | 2021 %+.1f\n",
            diff(ecB)[1], diff(ecB)[2], diff(dr)[1], diff(dr)[2]))
sink()
system(paste("cat", shQuote(file.path(OUT, "logs", "10_monthly_conceptions_both.txt"))))

saveRDS(list(ecA = ecA, ecB = ecB, dr = dr), file.path(OUT, "monthly_kl_both.rds"))

# ---- Two-panel comparative table: all conception years, percent only ---------
# (MM 2026-07-28: percent-only columns Mar-Dec 2020, 2021, 2022, 2023 (+2024 DR);
#  bold poorest-other differential row; counts version goes to the appendix.)
WINS <- list(
  "2020" = c(as.Date("2020-03-01"), as.Date("2020-12-01")),
  "2021" = c(as.Date("2021-01-01"), as.Date("2021-12-01")),
  "2022" = c(as.Date("2022-01-01"), as.Date("2022-12-01")),
  "2023" = c(as.Date("2023-01-01"), as.Date("2023-12-01")),
  "2024" = c(as.Date("2024-01-01"), as.Date("2024-12-01")))
pcts <- function(d, wins) vapply(wins, function(w) win(d, w[1], w[2])["pct"], numeric(1))
cells <- function(v, ncol) {
  out <- sprintf("%+.1f\\%%", v)
  if (length(out) < ncol) out <- c(out, rep("--", ncol - length(out)))
  paste(out, collapse = " & ")
}
row_line <- function(lab, d, wins, ncol, bold = FALSE) {
  s <- sprintf("%s & %s \\\\", lab, cells(pcts(d, wins), ncol))
  if (bold) s <- gsub("([+-][0-9.]+\\\\%)", "\\\\textbf{\\1}", s)
  s
}
dif_series <- function(lst, wins) pcts(lst$teen_poor, wins) - pcts(lst$teen_rich, wins)
panel2 <- function(title, lst, wins, ncol) c(
  sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", ncol + 1, title),
  row_line("All ages",            lst$all_ages,  wins, ncol),
  row_line("Teens 15--19 (all)",  lst$teen_all,  wins, ncol),
  row_line("\\quad Poorest 20\\%", lst$teen_poor, wins, ncol),
  row_line("\\quad Other 80\\%",   lst$teen_rich, wins, ncol),
  row_line("Women 20--29",        lst$adult2029, wins, ncol),
  sprintf("\\textbf{Poorest $-$ other, teens (pp)} & %s \\\\",
          {v <- dif_series(lst, wins); s <- sprintf("\\textbf{%+.1f}", v)
           if (length(s) < ncol) s <- c(s, "--"); paste(s, collapse = " & ")}),
  "\\addlinespace")
NC <- 5
tab <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Conceptions relative to pre-pandemic seasonal baselines, percent deviations by conception year}",
  "\\label{tab:conceptions_both}",
  sprintf("\\begin{tabular}{l%s}", paste(rep("c", NC), collapse = "")),
  "\\toprule",
  " & Mar--Dec 2020 & 2021 & 2022 & 2023 & 2024 \\\\",
  "\\midrule",
  panel2("Panel A. Ecuador (baseline: 2016--2019 seasonality and trend)", ecB, WINS[1:4], NC),
  panel2("Panel B. Dominican Republic (baseline: 2019 seasonal means)", dr, WINS, NC),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.95\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Conceptions dated as birth month minus nine months; predicted",
  "conceptions estimated out of sample from pre-pandemic data only, following",
  "Kearney and Levine (2023). Dominican baseline: month-of-year means on",
  "2019m1--2020m2 (the registration phase-in precludes earlier years; no trend",
  "is removed, so Dominican level deviations beyond 2021 increasingly reflect",
  "the secular decline in adolescent fertility, and the poorest-minus-other",
  "differential, which nets out components common to both groups, is the",
  "informative row). Ecuadorian baseline: month fixed effects plus linear",
  "trend on 2016m1--2020m2; the full 2010--2019 pre-period gives a 2020",
  "all-ages deviation of $-$7.4 percent instead of $-$9.4 percent. Ecuadorian",
  "conceptions are fully observed through 2023. Poorest 20\\% defined by the",
  "2010 census asset indices. Implied conception counts are reported in",
  "Appendix Table A1.",
  "\\end{minipage}",
  "\\end{table}")
writeLines(tab, file.path(OUT_TAB, "tab2_conceptions_both.tex"))

# Appendix counts table (2020/2021 windows with counts, as before)
row_tex <- function(lab, d) {
  b20 <- win(d, WINS[["2020"]][1], WINS[["2020"]][2]); y21 <- win(d, WINS[["2021"]][1], WINS[["2021"]][2])
  sprintf("%s & %+.1f\\%% & %s & %+.1f\\%% & %s \\\\", lab,
          b20["pct"], format(round(b20["n"]), big.mark = ","),
          y21["pct"], format(round(y21["n"]), big.mark = ","))
}
panelc <- function(title, lst) c(
  sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", title),
  row_tex("All ages",            lst$all_ages),
  row_tex("Teens 15--19 (all)",  lst$teen_all),
  row_tex("\\quad Poorest 20\\%", lst$teen_poor),
  row_tex("\\quad Other 80\\%",   lst$teen_rich),
  row_tex("Women 20--29",        lst$adult2029),
  "\\addlinespace")
taba <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Conception deviations with implied counts, 2020--2021}",
  "\\label{tab:conceptions_counts}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Mar--Dec 2020} & \\multicolumn{2}{c}{2021} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  " & \\% deviation & count & \\% deviation & count \\\\",
  "\\midrule",
  panelc("Panel A. Ecuador", ecB),
  panelc("Panel B. Dominican Republic", dr),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.94\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Conceptions are dated as birth month minus nine months. Deviations and",
  "implied counts compare observed conceptions with out-of-sample predictions from",
  "pre-pandemic seasonality and trends. The Dominican baseline uses month-of-year",
  "means estimated over January 2019--February 2020; the Ecuador baseline uses",
  "month fixed effects and a linear trend estimated over 2016--2019. Positive",
  "counts indicate excess conceptions and negative counts indicate shortfalls.",
  "Poorest 20\\% denotes the bottom quintile of the 2010 census asset index.",
  "\\end{minipage}",
  "\\end{table}")
writeLines(taba, file.path(OUT_TAB, "tabA1_conceptions_counts.tex"))

cat("\n== ALL-YEARS percent deviations ==\n")
show <- function(nm, lst, wins) for (g in names(lst))
  cat(sprintf("%s %-10s %s\n", nm, g, paste(sprintf("%+.1f", pcts(lst[[g]], wins)), collapse = "  ")))
show("DR", dr, WINS); show("EC", ecB, WINS[1:4])
cat("DR dif:", sprintf("%+.1f", dif_series(dr, WINS)), "\n")
cat("EC dif:", sprintf("%+.1f", dif_series(ecB, WINS[1:4])), "\n")
cat("Saved tab2_conceptions_both.tex and tabA1_conceptions_counts.tex\n")
