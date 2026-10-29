# 12_sensitivity.R — sensitivity analyses requested by OG (2026-07-29 meeting):
# (A) DR estimates using only registration-complete pre years (sample 2019-2025,
#     pre = 2019-2020, base 2020) -> Appendix Table A5.
# (B) Minority-concentration DDD: high-Haitian municipalities in the DR and
#     high-Indigenous cantons in Ecuador, defined as the PRE-pandemic (<=2019)
#     mean minority-mother share above the cross-locality mean (fixed, to avoid
#     endogeneity: the DR Haitian share doubles post-2020) -> Appendix Table A6.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

dr <- as.data.table(readRDS(file.path(OUT, "dr_panel.rds")))
ec <- as.data.table(haven::read_dta(file.path(EC_DATA, "base_canton_ENV_FE_rev2026.dta"),
        col_select = c("canton", "year", paste0("log_dif", 1:4), "post", "poorest",
                       "avgind", "avgage", "avgeduc", "sr", "avginst", "avgsba",
                       "avgcp", "avgchild2")))
ec <- ec[, lapply(.SD, as.numeric)]
XR_DR <- "avghaitian + avgage + avgeduc + sr + avginst + avgsba + avgcp"
XR_EC <- "avgind + avgage + avgeduc + sr + avginst + avgsba + avgcp + avgchild2"
stars <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .1, "*", "")))
pp <- function(m, term) coeftable(m)[term, ]

# ---- (A) DR, registration-complete sample 2019-2025 --------------------------
dr19 <- dr[year >= 2019]
a5 <- lapply(1:4, function(k) feols(
  as.formula(sprintf("log_dif%d ~ post:poorest + factor(year) + %s | id", k, XR_DR)),
  data = dr19, cluster = ~id))
cat("== A5: DR pooled DDD, sample 2019-2025 (vs full-sample estimate) ==\n")
drm <- readRDS(file.path(OUT, "dr_models_main.rds"))
rows5 <- character(0)
for (k in 1:4) {
  a <- pp(a5[[k]], "post:poorest")
  f <- pp(drm[[paste0("did_logdif", k, "_X")]], "post:poorest")
  cat(sprintf("vs %s: 2019+ %+.3f (%.3f)%s | full %+.3f (%.3f)%s\n", AG_LABS[AG_COMP[k]],
      a[1], a[2], stars(a[4]), f[1], f[2], stars(f[4])))
  rows5 <- c(rows5, sprintf("Women %s & %+.3f (%.3f)%s & %+.3f (%.3f)%s \\\\",
      AG_LABS[AG_COMP[k]], f[1], f[2], stars(f[4]), a[1], a[2], stars(a[4])))
}
writeLines(c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Dominican Republic: sensitivity to the registration-complete sample}",
  "\\label{tab:dr_2019pre}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Teen rate vs. & Full sample (2016--2025) & Registration-complete (2019--2025) \\\\",
  "\\midrule", rows5, "\\bottomrule", "\\end{tabular}",
  "\\begin{minipage}{0.85\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Post $\\times$ Poorest coefficients from the pooled DDD (log teen",
  "minus log comparison-group rate), municipality and year fixed effects,",
  "covariates, SEs clustered by municipality. The restricted sample drops the",
  "2016--2018 registration phase-in years, leaving pre = 2019--2020.",
  "*** p$<$.01, ** p$<$.05, * p$<$.1.",
  "\\end{minipage}", "\\end{table}"),
  file.path(OUT_TAB, "tabA5_dr_2019pre.tex"))

# ---- (B) Minority-concentration DDD ------------------------------------------
# Fixed pre-pandemic definition: mean minority-mother share over years <= 2019,
# above the cross-locality mean of those means.
mh <- dr[year <= 2019, .(hti_pre = mean(avghaitian, na.rm = TRUE)), by = id]
mh[, hi_min := as.integer(hti_pre > mean(hti_pre))]
dr <- merge(dr, mh[, .(id, hi_min)], by = "id")
mi <- ec[year <= 2019, .(ind_pre = mean(avgind, na.rm = TRUE)), by = canton]
mi[, hi_min := as.integer(ind_pre > mean(ind_pre))]
ec <- merge(ec, mi[, .(canton, hi_min)], by = "canton")
cat(sprintf("\nHigh-minority localities: DR %d/155 (Haitian), EC %d/221 (Indigenous)\n",
            uniqueN(dr[hi_min == 1, id]), uniqueN(ec[hi_min == 1, canton])))
cat(sprintf("Overlap with poorest quintile: DR %d/31, EC %d/44\n",
            uniqueN(dr[hi_min == 1 & poorest == 1, id]), uniqueN(ec[hi_min == 1 & poorest == 1, canton])))

a6 <- list()
rows6 <- character(0)
cat("\n== A6: minority-concentration pooled DDD (post x hi_minority) ==\n")
for (k in 1:4) {
  a6[[paste0("dr", k)]] <- feols(
    as.formula(sprintf("log_dif%d ~ post:hi_min + factor(year) + %s | id", k, XR_DR)),
    data = dr, cluster = ~id)
  a6[[paste0("ec", k)]] <- feols(
    as.formula(sprintf("log_dif%d ~ post:hi_min + factor(year) + %s | canton", k, XR_EC)),
    data = ec, cluster = ~canton)
  a <- pp(a6[[paste0("dr", k)]], "post:hi_min"); b <- pp(a6[[paste0("ec", k)]], "post:hi_min")
  cat(sprintf("vs %s: DR %+.3f (%.3f)%s | EC %+.3f (%.3f)%s\n", AG_LABS[AG_COMP[k]],
      a[1], a[2], stars(a[4]), b[1], b[2], stars(b[4])))
  rows6 <- c(rows6, sprintf("Women %s & %+.3f (%.3f)%s & %+.3f (%.3f)%s \\\\",
      AG_LABS[AG_COMP[k]], b[1], b[2], stars(b[4]), a[1], a[2], stars(a[4])))
}
writeLines(c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Minority-concentration sensitivity: Post $\\times$ high minority share}",
  "\\label{tab:minority}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Teen rate vs. & Ecuador: high Indigenous share & DR: high Haitian share \\\\",
  "\\midrule", rows6, "\\bottomrule", "\\end{tabular}",
  "\\begin{minipage}{0.9\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Pooled DDD (log teen minus log comparison-group rate) with Post",
  "$\\times$ high-minority-concentration dummy in place of the poorest-quintile",
  "dummy. High concentration = pre-pandemic ($\\leq$2019) mean share of births to",
  "Haitian mothers (DR) or Indigenous mothers (Ecuador) above the",
  "cross-locality mean. Locality and year fixed effects, covariates, clustered",
  "SEs. *** p$<$.01, ** p$<$.05, * p$<$.1.",
  "\\end{minipage}", "\\end{table}"),
  file.path(OUT_TAB, "tabA6_minority.tex"))
saveRDS(list(a5 = a5, a6 = a6), file.path(OUT, "sensitivity_models.rds"))
cat("\nSaved tabA5_dr_2019pre.tex and tabA6_minority.tex\n")
