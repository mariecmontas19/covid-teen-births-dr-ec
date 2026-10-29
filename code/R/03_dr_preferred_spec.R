# 03_dr_preferred_spec.R — Omar's preferred specification on DR data
# (log teen rate with comparison-group log rate as control + year x poorest),
# exactly matching his email regressions, plus adult-vs-adult placebos.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

panel <- as.data.table(readRDS(file.path(OUT, "dr_panel.rds")))
XVARS <- c("avghaitian", "avgage", "avgeduc", "sr", "avginst", "avgsba", "avgcp")
XRHS  <- paste(XVARS, collapse = " + ")

mods <- list()
# Omar's final spec: xtreg log_tn1519f log_tn{comp}f poorest i.year##poorest ib2020.year $X, fe vce(cluster)
for (g in AG_COMP) {
  mods[[paste0("pref_", g)]] <- feols(
    as.formula(sprintf(
      "log_tn1519 ~ log_tn%s + i(year, poorest, ref = 2020) + factor(year) + %s | id", g, XRHS)),
    data = panel, cluster = ~id)
}
# Simple DiD version of the same (his post##poorest with rate control)
for (g in AG_COMP) {
  mods[[paste0("pref_did_", g)]] <- feols(
    as.formula(sprintf("log_tn1519 ~ log_tn%s + post:poorest + factor(year) + %s | id", g, XRHS)),
    data = panel, cluster = ~id)
}
# Placebos: adult-vs-adult level differences (his tn_dif5-7), event-study form with X
for (k in 5:7) {
  mods[[paste0("placebo_tndif", k)]] <- feols(
    as.formula(sprintf("tn_dif%d ~ i(year, poorest, ref = 2020) + factor(year) + %s | id", k, XRHS)),
    data = panel, cluster = ~id)
}

saveRDS(mods, file.path(OUT, "dr_models_preferred.rds"))

sink(file.path(OUT, "logs", "03_dr_preferred_spec.txt"))
cat("== Omar-preferred spec: log teen rate ~ comparison log rate + year x poorest + X, muni FE ==\n")
print(etable(mods$pref_2024, mods$pref_2529, mods$pref_3034, mods$pref_3549, se.below = TRUE))
cat("\n== Same with post x poorest instead of year interactions ==\n")
print(etable(mods$pref_did_2024, mods$pref_did_2529, mods$pref_did_3034, mods$pref_did_3549, se.below = TRUE))
cat("\n== Placebos: adult-vs-adult rate differences, year x poorest ==\n")
print(etable(mods$placebo_tndif5, mods$placebo_tndif6, mods$placebo_tndif7, se.below = TRUE))
sink()

for (g in AG_COMP) {
  cf <- coef(mods[[paste0("pref_", g)]])
  post_cf <- cf[grep("year::202[1-5]:poorest", names(cf))]
  cat(sprintf("pref (teen | control %s) year x poorest 2021-2025: %s\n",
              AG_LABS[g], paste(sprintf("%.3f", post_cf), collapse = " ")))
}
cat("Log: analysis/output/logs/03_dr_preferred_spec.txt\n")
