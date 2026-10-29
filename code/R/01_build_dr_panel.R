# 01_build_dr_panel.R — DR municipality-year analysis panel, mirroring
# Omar's Ecuador canton panel (ENVanalysis_rev2026v2.log).
# Outcomes: birth rates per 1,000 women by age group (rateA denominator),
# log(rate+1), teen-minus-comparison differences; poorest = bottom 20% of
# the 2022-census wealth_index; post = year > 2020.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

# ---- 1. Rates by age group (rateA) from the dissertation panel -------------
ap <- as.data.table(readRDS(file.path(DR_CLEAN, "panel_age_long_155.rds")))
ap <- ap[ag %in% c("15-19", "20-24", "25-29", "30-34")]
ap[, agc := chartr("-", "_", ag)]
rates <- dcast(ap, id + year ~ agc, value.var = c("n_births", "womenA", "rateA"))
setnames(rates,
         c("rateA_15_19", "rateA_20_24", "rateA_25_29", "rateA_30_34"),
         c("tn1519", "tn2024", "tn2529", "tn3034"))

# ---- 2. Composite 35-49 group (births and Series A population, 155 geography) ------
b <- as.data.table(readRDS(file.path(DR_CLEAN, "births_clean_2016_2025.rds")))
b <- b[!is.na(adm3_pcode) & !is.na(birth_year)]
fold_map <- c(DOM010905 = "DOM010901", DOM051703 = "DOM051701", DOM012510 = "DOM012501")
b[, adm3_pcode := fifelse(adm3_pcode %in% names(fold_map),
                          unname(fold_map[adm3_pcode]), adm3_pcode)]
b3549 <- b[age_mom >= 35 & age_mom <= 49,
           .(nb_3549 = .N), by = .(id = adm3_pcode, year = birth_year)]

pop <- as.data.table(readRDS(file.path(DR_CLEAN, "pop_municipio_2016_2025_A.rds")))
w3549 <- pop[sex == "F" & age_group %in% c("35-39", "40-44", "45-49"),
             .(women3549 = sum(pop)), by = .(id = adm3_pcode, year)]
sr_dt <- pop[, .(pop_tot = sum(pop), pop_f = sum(pop[sex == "F"])),
             by = .(id = adm3_pcode, year)]
sr_dt[, sr := (pop_tot - pop_f) / pop_f * 100]

rates <- merge(rates, b3549,  by = c("id", "year"), all.x = TRUE)
rates <- merge(rates, w3549,  by = c("id", "year"), all.x = TRUE)
rates <- merge(rates, sr_dt[, .(id, year, sr)], by = c("id", "year"), all.x = TRUE)
rates[is.na(nb_3549), nb_3549 := 0]
rates[, tn3549 := nb_3549 / women3549 * 1000]

# ---- 3. Omar-analog time-varying covariates from births microdata ----------
# Per age group and muni-year: avg mother age, avg education (1-4), single-mother
# rate, institutional-delivery rate, skilled-attendant rate, 4+ ANC rate,
# Haitian-mother share (DR analog of Ecuador's Indigenous share).
# No birth-order variable exists in the DR microdata -> no child2 analog.
b[, agegrp := fcase(age_mom >= 15 & age_mom <= 19, "1519",
                    age_mom >= 20 & age_mom <= 24, "2024",
                    age_mom >= 25 & age_mom <= 29, "2529",
                    age_mom >= 30 & age_mom <= 34, "3034",
                    age_mom >= 35 & age_mom <= 49, "3549")]
cov <- b[!is.na(agegrp),
         .(age    = mean(age_mom, na.rm = TRUE),
           educ   = mean(as.integer(educ), na.rm = TRUE),
           single = mean(1L - in_union, na.rm = TRUE),
           inst   = mean(facility != "Home", na.rm = TRUE),
           sba    = mean(doctor_attended, na.rm = TRUE),
           cp     = mean(anc4, na.rm = TRUE),
           hti    = mean(nationality == "Haitian", na.rm = TRUE)),
         by = .(id = adm3_pcode, year = birth_year, agegrp)]
covw <- dcast(cov, id + year ~ agegrp, value.var = c("age", "educ", "single", "inst", "sba", "cp", "hti"))

panel <- merge(rates, covw, by = c("id", "year"), all.x = TRUE)

# Row-mean covariates across the five age groups (Omar's avg* variables),
# then impute remaining missings with municipality means (his egen/replace step)
grp5 <- c("1519", "2024", "2529", "3034", "3549")
for (v in c("age", "educ", "single", "inst", "sba", "cp", "hti")) {
  cols <- paste0(v, "_", grp5)
  panel[, (paste0("avg", v)) := rowMeans(.SD, na.rm = TRUE), .SDcols = cols]
}
setnames(panel, "avghti", "avghaitian")
XVARS <- c("avghaitian", "avgage", "avgeduc", "sr", "avginst", "avgsba", "avgcp")
for (v in XVARS) {
  panel[, (v) := fifelse(is.finite(get(v)), get(v), NA_real_)]
  panel[, (v) := fifelse(is.na(get(v)), mean(get(v), na.rm = TRUE), get(v)), by = id]
}

# ---- 4. Log outcomes and difference outcomes (Omar's construction) ---------
tn_all <- paste0("tn", c(AG_TEEN, AG_COMP))
panel[is.na(tn1519), tn1519 := 0]                      # his zero-fill of teen rate
for (v in tn_all) panel[, (paste0("log_", v)) := log(get(v) + 1)]
for (k in seq_along(AG_COMP)) {
  panel[, (paste0("tn_dif",  k)) := tn1519     - get(paste0("tn",     AG_COMP[k]))]
  panel[, (paste0("log_dif", k)) := log_tn1519 - get(paste0("log_tn", AG_COMP[k]))]
}
# Placebo differences among adult groups (his tn_dif5-7)
panel[, tn_dif5 := tn2024 - tn2529]
panel[, tn_dif6 := tn2024 - tn3034]
panel[, tn_dif7 := tn2024 - tn3549]

# ---- 5. Poorest quintile from the census wealth indices --------------------
# MAIN: 2010-census index (pre-pandemic, same vintage as Omar's Ecuador CPV2010
# index; built in 00b_wealth2010.R, KMO .926). Robustness: 2022-census index.
cc <- as.data.table(readRDS(file.path(DR_CLEAN, "census_covariates_muni.rds")))
cc <- cc[, .(id = adm3_pcode, adm3_name, prov_norm, wealth_index_2022 = wealth_index)]
w10 <- as.data.table(readRDS(file.path(OUT, "wealth2010_muni.rds")))
panel <- merge(panel, cc, by = "id", all.x = TRUE)
panel <- merge(panel, w10[, .(id = adm3_pcode, wealth_index_2010)], by = "id", all.x = TRUE)
u <- unique(panel[, .(id, wealth_index_2010, wealth_index_2022)])
panel[, poorest      := as.integer(wealth_index_2010 <= quantile(u$wealth_index_2010, .20, na.rm = TRUE))]
panel[, poorest_2022 := as.integer(wealth_index_2022 <= quantile(u$wealth_index_2022, .20, na.rm = TRUE))]
panel[, post := as.integer(year > 2020)]

stopifnot(uniqueN(panel$id) == 155, nrow(panel) == 155 * 10)
cat(sprintf("Panel: %d munis x %d years | poorest munis: %d | teen rate mean %.1f\n",
            uniqueN(panel$id), uniqueN(panel$year),
            uniqueN(panel[poorest == 1, id]), mean(panel$tn1519, na.rm = TRUE)))

saveRDS(panel, file.path(OUT, "dr_panel.rds"))
pd <- copy(panel)
setnames(pd, names(pd), gsub("[^A-Za-z0-9_]", "_", names(pd)))
haven::write_dta(pd, file.path(OUT, "dr_panel.dta"))
cat("Saved:", file.path(OUT, "dr_panel.rds"), "and .dta\n")
