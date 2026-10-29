# 14_mechanisms.R — mechanism regressions supported by the data audit
# (2026-07-29): (A) Ecuador mortality shock — canton-year deaths 2011-2022
# (EDG registry), log death rate on shock(2020-21) x poorest, canton + year FE.
# (B) DR infant mortality (2016-2023, unbalanced) as the DR-side suggestive
# check. Education and service channels are deferred: EC educ needs a code
# collapse around the 2018 instrument break; DR behavioral variables have no
# pre-period (educ from 2020, ANC from 2021, family planning from 2022).
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))
stars <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .1, "*", "")))

# ---- (A) Ecuador mortality -----------------------------------------------------
ed <- as.data.table(haven::read_dta(file.path(EC_DATA, "EDG", "EDG_canton_yr.dta"),
                                    col_select = c("year", "cant", "deaths")))
ed <- ed[, lapply(.SD, as.numeric)]
pop <- as.data.table(haven::read_dta(file.path(EC_DATA, "base_canton_ENV_FE_rev2026.dta"),
                                     col_select = c("canton", "year", "pob", "poorest", "hi_ind")))
pop <- pop[, lapply(.SD, as.numeric)]
ec <- merge(ed, pop, by.x = c("cant", "year"), by.y = c("canton", "year"))
cat(sprintf("EC mortality panel: %d canton-years, %d cantons, years %d-%d (dropped in merge: %d)\n",
            nrow(ec), uniqueN(ec$cant), min(ec$year), max(ec$year), nrow(ed) - nrow(ec)))
ec[, dr_all := deaths / pob * 1000]
ec[, shock := as.integer(year %in% 2020:2021)]
m_ec <- feols(log(dr_all) ~ shock:poorest | cant + year, data = ec, cluster = ~cant)
r <- coeftable(m_ec)["shock:poorest", ]
cat(sprintf("EC: shock(2020-21) x poorest on log death rate: %+.4f (%.4f)%s p=%.4f\n",
            r[1], r[2], stars(r[4]), r[4]))
# Also by year for the dynamics
m_ec_es <- feols(log(dr_all) ~ i(year, poorest, ref = 2019) | cant + year, data = ec, cluster = ~cant)
ct <- coeftable(m_ec_es)
ct20 <- ct[grep("year::202[0-2]:poorest", rownames(ct)), , drop = FALSE]
cat("EC year x poorest (ref 2019):", sprintf("%d %+.3f%s", 2020:2022, ct20[,1], stars(ct20[,4])), "\n")
# Indigenous-concentration margin (added 2026-07-30, OG's comment 5: does the
# mortality shock load on hi_ind rather than poorest, given their overlap?)
m_ec_ind <- feols(log(dr_all) ~ shock:hi_ind | cant + year, data = ec, cluster = ~cant)
ri <- coeftable(m_ec_ind)["shock:hi_ind", ]
cat(sprintf("EC: shock(2020-21) x hi_ind  on log death rate: %+.4f (%.4f)%s p=%.4f\n",
            ri[1], ri[2], stars(ri[4]), ri[4]))
m_ec_both <- feols(log(dr_all) ~ shock:poorest + shock:hi_ind | cant + year, data = ec, cluster = ~cant)
cb <- coeftable(m_ec_both)
cat(sprintf("EC joint: shock x poorest %+.4f (%.4f)%s | shock x hi_ind %+.4f (%.4f)%s\n",
            cb["shock:poorest",1], cb["shock:poorest",2], stars(cb["shock:poorest",4]),
            cb["shock:hi_ind",1],  cb["shock:hi_ind",2],  stars(cb["shock:hi_ind",4])))

# ---- (B) DR infant mortality ---------------------------------------------------
im <- as.data.table(readRDS(file.path(DR_CLEAN, "infant_mortality_muni_year.rds")))
w10 <- as.data.table(readRDS(file.path(OUT, "wealth2010_muni.rds")))
w10[, poorest := as.integer(wealth_index_2010 <= quantile(wealth_index_2010, .20))]
im <- merge(im, w10[, .(adm3_pcode, poorest)], by = "adm3_pcode")
cat(sprintf("\nDR IMR panel: %d muni-years, %d munis, years %d-%d\n",
            nrow(im), uniqueN(im$adm3_pcode), min(im$year), max(im$year)))
im[, shock := as.integer(year %in% 2020:2021)]
m_dr <- feols(log(IMR + 1) ~ shock:poorest | adm3_pcode + year,
              data = im[live_births > 0], cluster = ~adm3_pcode)
r2 <- coeftable(m_dr)["shock:poorest", ]
cat(sprintf("DR: shock(2020-21) x poorest on log(IMR+1): %+.4f (%.4f)%s p=%.4f\n",
            r2[1], r2[2], stars(r2[4]), r2[4]))

saveRDS(list(ec_mort = m_ec, ec_mort_es = m_ec_es, ec_mort_ind = m_ec_ind,
             ec_mort_both = m_ec_both, dr_imr = m_dr),
        file.path(OUT, "mechanism_models.rds"))
sink(file.path(OUT, "logs", "14_mechanisms.txt"))
cat("== EC mortality: log all-cause death rate, canton FE + year FE, cluster canton ==\n")
print(summary(m_ec)); print(summary(m_ec_es))
cat("\n== EC mortality: Indigenous-concentration margin (single + joint) ==\n")
print(summary(m_ec_ind)); print(summary(m_ec_both))
cat("\n== DR infant mortality: log(IMR+1), muni FE + year FE ==\n")
print(summary(m_dr))
sink()
cat("Saved mechanism models + log.\n")
