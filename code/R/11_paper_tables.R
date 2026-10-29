# 11_paper_tables.R — paper tables under the final numbering scheme
# (MM 2026-07-28): Table 1 descriptives (DR/EC x pre/post), Table 3 pooled DDD
# both countries, Table 4 alternative specs (rate-as-control + continuous
# gradient), Appendix A2 event-study coefficients, A3 placebos, A4 descriptives
# by wealth group. Table 2 (conceptions) is written by 10_monthly_conceptions_both.R.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

dr <- as.data.table(readRDS(file.path(OUT, "dr_panel.rds")))
ec <- as.data.table(haven::read_dta(file.path(EC_DATA, "base_canton_ENV_FE_rev2026.dta"),
        col_select = c("canton", "year", paste0("tn", c("1519","2024","2529","3034","3549"), "f"),
                       paste0("log_dif", 1:4), paste0("log_tn", c("2024","2529","3034","3549"), "f"),
                       "log_tn1519f", "tn_dif5", "tn_dif6", "tn_dif7", "post", "poorest",
                       "avgind", "avgage", "avgeduc", "sr", "avginst", "avgsba", "avgcp", "avgchild2")))
ec <- ec[, lapply(.SD, as.numeric)]
ai <- as.data.table(haven::read_dta(file.path(EC_DATA, "CPV2010", "CPV2010_aindex_canton.dta")))
ec <- merge(ec, ai[, .(canton = as.numeric(canton), asset_index)], by = "canton")
ue <- unique(ec[, .(canton, asset_index)])
ec[, pov_z := -(asset_index - mean(ue$asset_index)) / sd(ue$asset_index)]
XR_DR <- "avghaitian + avgage + avgeduc + sr + avginst + avgsba + avgcp"
XR_EC <- "avgind + avgage + avgeduc + sr + avginst + avgsba + avgcp + avgchild2"

# ============ TABLE 1 + A4: descriptives ======================================
ROWS <- list(
  c("Teen birth rate, 15--19",      "tn1519",  "tn1519f"),
  c("Birth rate, 20--24",           "tn2024",  "tn2024f"),
  c("Birth rate, 25--29",           "tn2529",  "tn2529f"),
  c("Birth rate, 30--34",           "tn3034",  "tn3034f"),
  c("Birth rate, 35--49",           "tn3549",  "tn3549f"),
  c("Mother's mean age",            "avgage",  "avgage"),
  c("Mother's mean education$^a$",  "avgeduc", "avgeduc"),
  c("Institutional delivery",       "avginst", "avginst"),
  c("Skilled birth attendant",      "avgsba",  "avgsba"),
  c("4+ antenatal visits",          "avgcp",   "avgcp"),
  c("Sex ratio",                    "sr",      "sr"),
  c("Minority mother share$^b$",    "avghaitian", "avgind"))
msd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
msd2 <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
cell <- function(d, v, small) if (small) msd2(d[[v]]) else msd(d[[v]])
# DR pre column = 2019-2020 only (registration-complete years; MM 2026-07-29).
# Country order: Ecuador first everywhere (MM 2026-08-03).
desc_cols <- function(drd, ecd) {
  lines <- character(0)
  for (r in ROWS) {
    small <- r[2] %in% c("avginst", "avgsba", "avgcp", "avghaitian")
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s \\\\", r[1],
      cell(ecd[post == 0], r[3], small), cell(ecd[post == 1], r[3], small),
      cell(drd[post == 0 & year >= 2019], r[2], small), cell(drd[post == 1], r[2], small)))
  }
  lines
}
struct_rows <- function(drd, ecd) c(
  sprintf("Localities & \\multicolumn{2}{c}{%d} & \\multicolumn{2}{c}{%d} \\\\",
          uniqueN(ecd$canton), uniqueN(drd$id)),
  sprintf("Poorest-quintile localities & \\multicolumn{2}{c}{%d} & \\multicolumn{2}{c}{%d} \\\\",
          uniqueN(ecd[poorest == 1, canton]), uniqueN(drd[poorest == 1, id])),
  sprintf("Locality-years & %d & %d & %d & %d \\\\",
          nrow(ecd[post == 0]), nrow(ecd[post == 1]),
          nrow(drd[post == 0 & year >= 2019]), nrow(drd[post == 1])))
tab1 <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Descriptive statistics, before and after the pandemic}",
  "\\label{tab:descriptives}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Ecuador} & \\multicolumn{2}{c}{Dominican Republic} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  " & Pre (2010--20) & Post (2021--24) & Pre (2019--20) & Post (2021--25) \\\\",
  "\\midrule",
  desc_cols(dr, ec),
  "\\midrule",
  struct_rows(dr, ec),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.92\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Unweighted locality-year means; standard deviations in parentheses.",
  "Birth rates per 1,000 women in the age group. Post = birth years 2021 and",
  "later; birth-year 2020 is pre-pandemic because its births were conceived",
  "April 2019--March 2020 (only 4.43 percent of Dominican 2020 births were",
  "conceived after the March 19, 2020 state of emergency). The Dominican pre",
  "column uses 2019--2020 only, the registration-complete years; estimation",
  "samples use 2016--2025 with year fixed effects. $^a$Education coded on",
  "national scales (DR 1--4; Ecuador 0--8), not comparable across countries.",
  "$^b$Haitian-mother share in the DR; Indigenous-mother share in Ecuador.",
  "The DR national population-weighted teen rate peaks at 61.8 per 1,000 in",
  "2019. Descriptives by wealth quintile appear in Appendix Table A4.",
  "\\end{minipage}",
  "\\end{table}")
writeLines(tab1, file.path(OUT_TAB, "tab1_descriptives.tex"))

tabA4 <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Descriptive statistics by wealth group}",
  "\\label{tab:desc_wealth}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Ecuador} & \\multicolumn{2}{c}{Dominican Republic} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  " & Pre & Post & Pre & Post \\\\",
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Panel A. Poorest quintile}} \\\\",
  desc_cols(dr[poorest == 1], ec[poorest == 1]),
  "\\addlinespace",
  "\\multicolumn{5}{l}{\\textit{Panel B. Other localities}} \\\\",
  desc_cols(dr[poorest == 0], ec[poorest == 0]),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.96\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Entries are unweighted locality-year means, with standard deviations in",
  "parentheses. Birth rates are per 1,000 women in the indicated age group. The",
  "poorest quintile is defined using the 2010 census asset index. For Ecuador,",
  "Pre = 2010--2020 and Post = 2021--2024; for the Dominican Republic,",
  "Pre = 2019--2020 and Post = 2021--2025. $^{a}$Education uses national coding",
  "scales and is not directly comparable across countries. $^{b}$Minority share",
  "is the share of births to Haitian mothers in the Dominican Republic and to",
  "Indigenous mothers in Ecuador.",
  "\\end{minipage}",
  "\\end{table}")
writeLines(tabA4, file.path(OUT_TAB, "tabA4_descriptives_by_wealth.tex"))

# Markdown version of Table 1 for the Word draft
mdcell <- function(d, v, small) if (small) msd2(d[[v]]) else msd(d[[v]])
cat("\n== TABLE 1 (markdown) ==\n")
for (r in ROWS) {
  small <- r[2] %in% c("avginst", "avgsba", "avgcp", "avghaitian")
  cat(sprintf("| %s | %s | %s | %s | %s |\n", gsub("\\$.*?\\$", "*", r[1]),
    mdcell(dr[post == 0 & year >= 2019], r[2], small), mdcell(dr[post == 1], r[2], small),
    mdcell(ec[post == 0], r[3], small), mdcell(ec[post == 1], r[3], small)))
}
cat(sprintf("| Localities | %d | | %d | |\n| Poorest-quintile localities | %d | | %d | |\n| Locality-years | %d | %d | %d | %d |\n",
    uniqueN(dr$id), uniqueN(ec$canton), uniqueN(dr[poorest == 1, id]), uniqueN(ec[poorest == 1, canton]),
    nrow(dr[post == 0 & year >= 2019]), nrow(dr[post == 1]), nrow(ec[post == 0]), nrow(ec[post == 1])))

# ============ TABLE 3: pooled DDD both countries ==============================
drm <- readRDS(file.path(OUT, "dr_models_main.rds"))
ec_ddd <- lapply(1:4, function(k) feols(
  as.formula(sprintf("log_dif%d ~ post:poorest + factor(year) | canton", k)),
  data = ec, cluster = ~canton))
saveRDS(ec_ddd, file.path(OUT, "ec_models_ddd_pooled.rds"))
get_pp <- function(m, term = "post:poorest") coeftable(m)[term, ]
stars <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .1, "*", "")))
cat("\n== TABLE 3 rows (DDD) ==\n")
for (k in 1:4) {
  a <- get_pp(drm[[paste0("did_logdif", k)]]); b <- get_pp(ec_ddd[[k]])
  cat(sprintf("| Women %s | %+.3f (%.3f)%s | %+.3f (%.3f)%s |\n", AG_LABS[AG_COMP[k]],
      a[1], a[2], stars(a[4]), b[1], b[2], stars(b[4])))
}

# ============ TABLE 4: alternative specifications =============================
drp <- readRDS(file.path(OUT, "dr_models_preferred.rds"))
rob <- readRDS(file.path(OUT, "dr_models_robustness.rds"))
ec_rc <- lapply(AG_COMP, function(g) feols(
  as.formula(sprintf("log_tn1519f ~ log_tn%sf + post:poorest + factor(year) + %s | canton", g, XR_EC)),
  data = ec, cluster = ~canton))
saveRDS(ec_rc, file.path(OUT, "ec_models_ratecontrol_pooled.rds"))
cat("\n== TABLE 4 rows (rate-as-control | continuous gradient) ==\n")
t4 <- character(0)
for (k in 1:4) {
  a <- get_pp(drp[[paste0("pref_did_", AG_COMP[k])]])
  b <- get_pp(ec_rc[[k]])
  c_ <- get_pp(rob[[paste0("did_cont", k)]], "post:pov_z")
  d_ <- get_pp(rob[[paste0("ec_did_cont", k)]], "post:pov_z")
  line <- sprintf("| Women %s | %+.3f (%.3f)%s | %+.3f (%.3f)%s | %+.3f (%.3f)%s | %+.3f (%.3f)%s |",
      AG_LABS[AG_COMP[k]], b[1], b[2], stars(b[4]), a[1], a[2], stars(a[4]),
      d_[1], d_[2], stars(d_[4]), c_[1], c_[2], stars(c_[4]))
  cat(line, "\n"); t4 <- c(t4, gsub("^\\| | \\|$", "", gsub(" \\| ", " & ", line)))
}
tab4 <- c(
  "\\begin{table}[H]\\centering",  # [H]: MM wants main-text exhibits pinned in place (2026-07-31)
  "\\caption{Alternative specifications: comparison-group rate as control and continuous poverty gradient}",
  "\\label{tab:altspecs}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Rate as control: Post $\\times$ Poorest} & \\multicolumn{2}{c}{Continuous: Post $\\times$ Poverty (1 SD)} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Teen rate vs. & Ecuador & Dominican Rep. & Ecuador & Dominican Rep. \\\\",
  "\\midrule",
  paste0(t4, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.92\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Columns 1--2: log teen birth rate regressed on the comparison-group",
  "log rate, Post $\\times$ Poorest, locality and year fixed effects, and",
  "covariates. Columns 3--4: log difference (teen minus comparison group)",
  "regressed on Post $\\times$ standardized poverty index (higher = poorer).",
  "Standard errors clustered at the locality level. *** p$<$.01, ** p$<$.05, * p$<$.1.",
  "\\end{minipage}",
  "\\end{table}")
writeLines(tab4, file.path(OUT_TAB, "tab4_altspecs_both.tex"))

# ============ APPENDIX A2: event-study coefficients ===========================
ecm <- readRDS(file.path(OUT, "ec_models_poorest.rds"))
esrow <- function(m, yrs) {
  ct <- coeftable(m); out <- character(0)
  for (y in yrs) {
    r <- ct[grep(sprintf("year::%d:poorest", y), rownames(ct)), , drop = FALSE]
    out <- c(out, if (nrow(r)) sprintf("%+.3f (%.3f)%s", r[1,1], r[1,2], stars(r[1,4])) else "--")
  }
  out
}
a2 <- c("\\begin{table}[htbp]\\centering",
  "\\caption{Triple-difference event-study coefficients (year $\\times$ poorest)}",
  "\\label{tab:es_coefs}", "\\begin{tabular}{lcccc}", "\\toprule",
  "Year & vs. 20--24 & vs. 25--29 & vs. 30--34 & vs. 35--49 \\\\", "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Panel A. Ecuador (reference 2020)}} \\\\")
for (y in c(2010:2019, 2021:2024)) a2 <- c(a2, sprintf("%d & %s \\\\", y,
  paste(vapply(1:4, function(k) esrow(ecm[[paste0("es_logdif", k)]], y), character(1)), collapse = " & ")))
a2 <- c(a2, "\\addlinespace", "\\multicolumn{5}{l}{\\textit{Panel B. Dominican Republic (reference 2020)}} \\\\")
for (y in c(2016:2019, 2021:2025)) a2 <- c(a2, sprintf("%d & %s \\\\", y,
  paste(vapply(1:4, function(k) esrow(drm[[paste0("es_logdif", k)]], y), character(1)), collapse = " & ")))
a2 <- c(a2, "\\bottomrule", "\\end{tabular}",
  "\\begin{minipage}{0.94\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Each column reports coefficients from a separate triple-difference",
  "event-study regression. The outcome is the log teen (15--19) birth rate minus",
  "the log birth rate of the indicated comparison group. Entries are year",
  "$\\times$ poorest-quintile coefficients relative to 2020. All models include",
  "locality and year fixed effects without registry-derived covariates. Standard errors,",
  "clustered by municipality in the Dominican Republic and canton in Ecuador, are",
  "in parentheses. *** p$<$.01, ** p$<$.05, * p$<$.10.",
  "\\end{minipage}", "\\end{table}")
writeLines(a2, file.path(OUT_TAB, "tabA2_eventstudy_coefs.tex"))

# ============ APPENDIX A3: placebos (pooled) ==================================
dr[, `:=`(adult_logdif_2529 = log_tn2024 - log_tn2529,
          adult_logdif_3034 = log_tn2024 - log_tn3034,
          adult_logdif_3549 = log_tn2024 - log_tn3549)]
ec[, `:=`(adult_logdif_2529 = log_tn2024f - log_tn2529f,
          adult_logdif_3034 = log_tn2024f - log_tn3034f,
          adult_logdif_3549 = log_tn2024f - log_tn3549f)]
plc <- function(dt, outcome, id) feols(
  as.formula(sprintf("%s ~ post:poorest + factor(year) | %s", outcome, id)),
  data = dt, cluster = as.formula(paste0("~", id)))
a3rows <- character(0)
comps <- c("2529", "3034", "3549")
labs <- c("2529" = "20--24 vs.\\ 25--29", "3034" = "20--24 vs.\\ 30--34", "3549" = "20--24 vs.\\ 35--49")
cat("\n== A3 placebos (post x poorest, log differences) ==\n")
for (g in comps) {
  y <- paste0("adult_logdif_", g)
  a <- get_pp(plc(dr, y, "id")); b <- get_pp(plc(ec, y, "canton"))
  cat(sprintf("%s: DR %+.3f (%.3f) p=%.3f | EC %+.3f (%.3f) p=%.3f\n", g, a[1], a[2], a[4], b[1], b[2], b[4]))
  a3rows <- c(a3rows, sprintf("%s & %+.3f (%.3f)%s & %+.3f (%.3f)%s \\\\",
      labs[g], b[1], b[2], stars(b[4]), a[1], a[2], stars(a[4])))
}
a3 <- c("\\begin{table}[htbp]\\centering",
  "\\caption{Placebo contrasts among adult age groups (Post $\\times$ Poorest)}",
  "\\label{tab:placebos}", "\\begin{tabular}{lcc}", "\\toprule",
  "Log-rate difference & Ecuador & Dominican Republic \\\\", "\\midrule", a3rows,
  "\\bottomrule", "\\end{tabular}",
  "\\begin{minipage}{0.8\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Outcome is the log(rate+1) difference between ages 20--24 and the",
  "indicated older group. Models include locality and unrestricted year fixed",
  "effects without registry-derived covariates; SEs are clustered by locality.",
  "\\end{minipage}", "\\end{table}")
writeLines(a3, file.path(OUT_TAB, "tabA3_placebos.tex"))

# ============ APPENDIX A7: rate-as-control event-study coefficients ===========
# OG's specification (log teen rate ~ comparison log rate + i.year##poorest + X),
# year by year. The Ecuador columns reproduce his ENVanalysis_rev2026v2.log
# poorest regressions verbatim (verified to 4 decimals, 2026-07-29); his log has
# no 30-34 regression, so that column is our extension of his spec.
drp <- readRDS(file.path(OUT, "dr_models_preferred.rds"))
GRP <- c("2024", "2529", "3034", "3549")
a7 <- c("\\begin{table}[htbp]\\centering",
  "\\caption{Event-study coefficients from the rate-as-control specification (year $\\times$ poorest)}",
  "\\label{tab:es_ratecontrol}", "\\begin{tabular}{lcccc}", "\\toprule",
  "Year & vs. 20--24 & vs. 25--29 & vs. 30--34 & vs. 35--49 \\\\", "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Panel A. Ecuador (reference 2020)}} \\\\")
for (y in c(2010:2019, 2021:2024)) a7 <- c(a7, sprintf("%d & %s \\\\", y,
  paste(vapply(GRP, function(g) esrow(ecm[[paste0("pref_", g)]], y), character(1)), collapse = " & ")))
a7 <- c(a7, "\\addlinespace", "\\multicolumn{5}{l}{\\textit{Panel B. Dominican Republic (reference 2020)}} \\\\")
for (y in c(2016:2019, 2021:2025)) a7 <- c(a7, sprintf("%d & %s \\\\", y,
  paste(vapply(GRP, function(g) esrow(drp[[paste0("pref_", g)]], y), character(1)), collapse = " & ")))
a7 <- c(a7, "\\bottomrule", "\\end{tabular}",
  "\\begin{minipage}{0.94\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Alternative specification to Appendix Table~\\ref{tab:es_coefs}. The",
  "outcome is the log teen (15--19) birth rate, regressed on the log birth rate",
  "of the indicated comparison group (entering as a control with a freely",
  "estimated coefficient) rather than differencing it out. Entries are year",
  "$\\times$ poorest-quintile coefficients relative to 2020. All models include",
  "locality and year fixed effects and time-varying covariates. Standard errors,",
  "clustered by municipality in the Dominican Republic and canton in Ecuador, are",
  "in parentheses. *** p$<$.01, ** p$<$.05, * p$<$.10.",
  "\\end{minipage}", "\\end{table}")
writeLines(a7, file.path(OUT_TAB, "tabA7_ratecontrol_es.tex"))

# ============ APPENDIX A8: KMO sampling adequacy for the asset indices ========
# Panel A (Ecuador): values transcribed from the Stata construction log of the
# CPV2010 canton index, covid_baby_bust_boom/logs/cpv2010_gini_v2.log
# (`estat kmo` after `pca high_occup H07-H11 ownkitchen ownshower cleancooking
# drinkw1 rooft1 wall_wood floor_parquet_tile water_piped_inside toilet_sewer`,
# N = 3,815,527 households, overall KMO 0.9014). Item labels from
# violence/dofiles/ecuador_asset_gini_2010.do.
# Panel B (DR): computed in 00b_wealth2010.R (psych::KMO on the final retained
# item set) and saved to wealth2010_kmo.rds.
ec_kmo <- list(
  c("Persons per bedroom above national mean", "0.924"),
  c("Telephone land line",                     "0.936"),
  c("Cell phone",                              "0.900"),
  c("Internet service",                        "0.826"),
  c("Computer",                                "0.863"),
  c("Cable or satellite television",           "0.934"),
  c("Separate kitchen space",                  "0.910"),
  c("Separate shower space",                   "0.914"),
  c("Clean cooking fuel",                      "0.882"),
  c("Drinking water consumed as received",     "0.880"),
  c("Roof is cement",                          "0.919"),
  c("Walls are wood",                          "0.927"),
  c("Floor is parquet or tile",                "0.923"),
  c("Piped water inside the dwelling",         "0.899"),
  c("Toilet connected to public sewer",        "0.904"))
DR_LABS <- c(
  H09A_ESTUFA = "Stove", H09B_NEVERA = "Refrigerator", H09C_LAVADORA = "Washing machine",
  H09D_TELEVISOR = "Television", H09E_RADIO = "Radio", H09F_CISTERNA = "Water cistern",
  H09G_TINACO = "Rooftop water tank", H09H_COMPUTADORA = "Computer",
  H09I_INTERNET = "Internet service", H09J_INVERSOR = "Power inverter",
  H09K_PLANTA = "Electric generator", H09L_TELEFONO_FIJO = "Telephone land line",
  H09M_CELULAR = "Cell phone", H09N_AIRE = "Air conditioning",
  H09NN_MOTOR = "Motorcycle", H09O_AUTOMOVIL = "Automobile",
  crowded = "More than two persons per bedroom", toilet_flush = "Flush toilet",
  toilet_sole = "Toilet not shared with other households",
  garbage_coll = "Garbage collected by municipal service",
  water_piped_in = "Piped water inside the dwelling",
  cook_gas_elec = "Cooks with propane gas or electricity",
  light_grid = "Electric lighting from public grid")
kdr <- readRDS(file.path(OUT, "wealth2010_kmo.rds"))
lab_dr <- ifelse(is.na(DR_LABS[kdr$items]), kdr$items, DR_LABS[kdr$items])
a8 <- c("\\begin{table}[htbp]\\centering",
  "\\caption{Kaiser--Meyer--Olkin measures of sampling adequacy, 2010-census asset indices}",
  "\\label{tab:kmo}", "\\begin{tabular}{lc}", "\\toprule",
  "Item & KMO \\\\", "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Panel A. Ecuador (CPV 2010)}} \\\\",
  vapply(ec_kmo, function(r) sprintf("%s & %s \\\\", r[1], r[2]), character(1)),
  "\\textbf{Overall} & \\textbf{0.901} \\\\",
  "\\addlinespace",
  "\\multicolumn{2}{l}{\\textit{Panel B. Dominican Republic (Censo 2010)}} \\\\",
  sprintf("%s & %.3f \\\\", lab_dr, kdr$msai),
  sprintf("\\textbf{Overall} & \\textbf{%.3f} \\\\", kdr$overall),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{minipage}{0.78\\textwidth}\\vspace{2pt}\\footnotesize",
  "Notes: Kaiser--Meyer--Olkin measures of sampling adequacy for the household",
  "asset items entering each country's principal-components index, computed on",
  "the 2010 census household microdata (Ecuador: 3,815,527 households;",
  sprintf("Dominican Republic: %s complete-case households).",
          format(kdr$n_complete, big.mark = ",")),
  "In the Dominican Republic, items with individual measures below 0.85 (stove,",
  "internet service, motorcycle) were dropped before estimation; the Ecuadorian",
  "index retains its original 15-item set.",
  "\\end{minipage}", "\\end{table}")
writeLines(a8, file.path(OUT_TAB, "tabA8_kmo.tex"))
cat("\nAll paper tables written to", OUT_TAB, "\n")
