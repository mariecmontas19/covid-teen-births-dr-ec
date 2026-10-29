# 04_ecuador_replicate.R — replicate Omar's final Ecuador regressions
# (ENVanalysis_rev2026v2.log, poorest-quintile specs) in R with fixest,
# from his saved analysis file base_canton_ENV_FE_rev2026.dta.
# Purpose: (i) validate that feols reproduces xtreg,fe vce(cluster) exactly,
# (ii) store the Ecuador models for the comparative figures/tables.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

ec <- as.data.table(haven::read_dta(file.path(EC_DATA, "base_canton_ENV_FE_rev2026.dta")))
XRHS <- "avgind + avgage + avgeduc + sr + avginst + avgsba + avgcp + avgchild2"

pref <- function(comp) feols(
  as.formula(sprintf(
    "log_tn1519f ~ log_tn%sf + i(year, poorest, ref = 2020) + factor(year) + %s | canton",
    comp, XRHS)),
  data = ec, cluster = ~canton)

mods <- list(pref_2024 = pref("2024"), pref_2529 = pref("2529"), pref_3549 = pref("3549"))
# The missing piece in his log: the 30-34 comparison (email item (c) labeling issue)
mods$pref_3034 <- pref("3034")
# Headline-format analogs for the comparative figures: log-difference event studies
for (k in 1:4) {
  mods[[paste0("es_logdif", k)]] <- feols(
    as.formula(sprintf("log_dif%d ~ i(year, poorest, ref = 2020) + factor(year) | canton", k)),
    data = ec, cluster = ~canton)
  mods[[paste0("es_logdif", k, "_X")]] <- feols(
    as.formula(sprintf("log_dif%d ~ i(year, poorest, ref = 2020) + factor(year) + %s | canton", k, XRHS)),
    data = ec, cluster = ~canton)
}
saveRDS(mods, file.path(OUT, "ec_models_poorest.rds"))

# --- Verification against the printed Stata log (post-2020 year x poorest) ---
target <- list(
  pref_2024 = c(`2021` = .1331238, `2022` = .1652499, `2023` = .1363308, `2024` = .1612521),
  pref_2529 = c(`2021` = .2009409, `2022` = .2158326, `2023` = .1955722, `2024` = .2114764),
  pref_3549 = c(`2021` = .1612304, `2022` = .1876374, `2023` = .1707622, `2024` = .1708332))
ok <- TRUE
for (m in names(target)) {
  cf <- coef(mods[[m]])
  got <- cf[sprintf("year::%s:poorest", names(target[[m]]))]
  diff <- max(abs(got - target[[m]]))
  cat(sprintf("%s: max |R - Stata| = %.2e %s\n", m, diff, if (diff < 5e-6) "OK" else "MISMATCH"))
  if (diff >= 5e-6) ok <- FALSE
}
cf34 <- coef(mods$pref_3034)
cat("NEW pref_3034 (missing from his log) year x poorest 2021-2024:",
    sprintf("%.3f", cf34[grep("year::202[1-4]:poorest", names(cf34))]), "\n")
stopifnot(ok)
cat("Ecuador replication verified: fixest matches xtreg to <5e-6.\n")
