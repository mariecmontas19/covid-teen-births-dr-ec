# 15_household_margin.R — household-margin analyses approved by MM (2026-07-29):
# (1) Subsidized-insurance DDD: the mother's SeNaSa subsidized regime is a
#     means-tested household-poverty proxy attached to each DR birth from 2018.
#     Test: did teen births shift toward subsidized (poor) households after the
#     pandemic, relative to adult births, within municipalities? Cells are
#     municipality x year x teen x subsidized; log-count triple difference plus
#     a composition event study; split by poorest vs other municipalities.
# (2) Province-level concentration index: demonstration that coarser geography
#     attenuates the DR spatial gradient (footnote for the paper).
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))
stars <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .1, "*", "")))

b <- as.data.table(readRDS(file.path(DR_CLEAN, "births_clean_2016_2025.rds")))
b <- b[birth_year >= 2018 & !is.na(adm3_pcode) & !is.na(subsidized) & !is.na(age_mom)]
b[, teen := as.integer(age_mom >= 15 & age_mom <= 19)]
b <- b[teen == 1 | (age_mom >= 20 & age_mom <= 34)]
cat(sprintf("Analysis births 2018-2025 with insurance info: %s (%.1f%% of 15-34 births)\n",
            format(nrow(b), big.mark = ","), 100 * nrow(b) /
            nrow(readRDS(file.path(DR_CLEAN, "births_clean_2016_2025.rds")))))
w10 <- as.data.table(readRDS(file.path(OUT, "wealth2010_muni.rds")))
w10[, poorest := as.integer(wealth_index_2010 <= quantile(wealth_index_2010, .20))]
b <- merge(b, w10[, .(adm3_pcode, poorest)], by = "adm3_pcode")

# ---- (1a) log-count triple difference ---------------------------------------
cells <- b[, .N, by = .(id = adm3_pcode, year = birth_year, teen, sub = subsidized, poorest)]
cells <- cells[CJ(id = unique(id), year = 2018:2025, teen = 0:1, sub = 0:1), on = .(id, year, teen, sub)]
cells[is.na(N), N := 0]
cells[, poorest := poorest[!is.na(poorest)][1], by = id]
cells[, post := as.integer(year > 2020)]
run3 <- function(d, lab) {
  m <- feols(log(N + 1) ~ teen * sub * post | id + year, data = d, cluster = ~id)
  r <- coeftable(m)["teen:sub:post", ]
  cat(sprintf("%-22s teen x sub x post: %+.4f (%.4f)%s p=%.4f\n", lab, r[1], r[2], stars(r[4]), r[4]))
  m
}
cat("\n== (1a) Log-count DDD: did teen births shift toward subsidized households post-2020? ==\n")
m_all  <- run3(cells, "All municipalities")
m_poor <- run3(cells[poorest == 1], "Poorest quintile")
m_rich <- run3(cells[poorest == 0], "Other 80% (mixed/rich)")

# ---- (1b) composition event study: sub share of teen minus adult births ------
comp <- b[, .(sub_share = mean(subsidized)), by = .(id = adm3_pcode, year = birth_year, teen)]
comp <- dcast(comp, id + year ~ teen, value.var = "sub_share")
setnames(comp, c("0", "1"), c("sub_adult", "sub_teen"))
comp <- comp[!is.na(sub_teen) & !is.na(sub_adult)]
comp[, D := sub_teen - sub_adult]
m_es <- feols(D ~ i(year, ref = 2019) | id, data = comp, cluster = ~id)
ct <- coeftable(m_es)
cat("\n== (1b) Teen-minus-adult subsidized share, event study (ref 2019) ==\n")
for (y in c(2018, 2020:2025)) {
  r <- ct[sprintf("year::%d", y), ]
  cat(sprintf("%d: %+.4f (%.4f)%s\n", y, r[1], r[2], stars(r[4])))
}
cat(sprintf("Baseline gap 2019 (teen minus adult subsidized share): %+.3f\n",
            comp[year == 2019, mean(D)]))
saveRDS(list(ddd_all = m_all, ddd_poor = m_poor, ddd_rich = m_rich, comp_es = m_es),
        file.path(OUT, "household_margin_models.rds"))

# ---- (2) province-level concentration index (attenuation demo) ---------------
dr <- as.data.table(readRDS(file.path(OUT, "dr_panel.rds")))
conc_index <- function(d) {
  d <- d[order(wealth)]
  w <- d$women / sum(d$women); R <- cumsum(w) - w / 2
  y <- d$births / d$women; mu <- sum(w * y)
  sum(w * y * (2 * R - 1)) / mu
}
agg <- function(dd, lvl) dd[, .(births = sum(n_births_15_19), women = sum(womenA_15_19),
                                wealth = weighted.mean(wealth_index_2010, womenA_15_19)),
                            by = lvl]
cat("\n== (2) DR concentration index by aggregation level ==\n")
for (per in list(c("Pre (2019-20)", 2019, 2020), c("Post (2021-25)", 2021, 2025))) {
  dd <- dr[year >= as.integer(per[2]) & year <= as.integer(per[3])]
  cat(sprintf("%-14s municipality C = %+.4f | province C = %+.4f\n", per[1],
              conc_index(agg(dd, "id")), conc_index(agg(dd, "prov_norm"))))
}
cat("Done.\n")
