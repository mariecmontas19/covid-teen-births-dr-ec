# 16_census_gradient.R — household-level teen-birth wealth gradient from the
# 2022 censuses (OG's EQ_teen_birth_gradient_07_2026.do design, corrected).
# ECUADOR: replicate OG's headline gradient from his merged HH file, then build
#   properly WINDOWED outcomes. Note on his filters (Stata missing = +infinity):
#   `if latest_teen_birth_yr>2020` keeps never-teen-birth HHs (intended as the
#   denominator) but drops pre-2021 teen-birth HHs; his "pre" filter
#   `<=2020 | >2015` is TRUE for all HHs, so his pre index is the all-history
#   gradient. Corrected windows here: pre = any teen first birth 2016-2020,
#   post = any in 2021-2022; denominators = all households.
#   Windowed pre indicator is exact for HHs with <=2 teen births (earliest/
#   latest logic); ambiguous >=3-birth straddlers are counted and reported.
# DR: analog from the 2022 census (BD_FINAL_VIVIENDA_HOGAR_PERSONA):
#   HH wealth PCA (dissertation v2 recipe), quintiles; teen birth dated from
#   the LAST live birth (P67) and mother's age (DR census lacks age at first
#   birth), so pre-window counts are eroded where a later birth overwrote a
#   2016-20 teen birth — stated caveat; post window (2021-22) is clean.
# Output: fig5 (X = Poorest..Richest, Y = % of HH per year, pre vs post),
#   gradient slopes, Erreygers-corrected concentration indices.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

QLAB <- c("Poorest", "Poorer", "Middle", "Richer", "Richest")
err_ci <- function(y, rank_var) {          # Erreygers normalized CI (binary y)
  o <- order(rank_var)
  y <- y[o]; n <- length(y)
  R <- (seq_len(n) - 0.5) / n
  mu <- mean(y)
  C <- 2 * mean(y * R) / mu - 1
  4 * mu * C
}

sink(file.path(OUT, "logs", "16_census_gradient.txt"), split = TRUE)

# ============================ ECUADOR ========================================
# Rebuild the full-sample HH file from OG's components (his 1:1 ID_HOG merge).
# Do NOT read cpv2022_hh_all.dta: his do-file re-saves it at the END of the
# script, AFTER `keep if latest_teen_birth_yr>2015`, so its content depends on
# where his last run stopped (2026-07-30 sync = filtered 3.8M-HH snapshot;
# his logged regression N = 5,188,827 uses the pre-keep full sample).
# Merge key must be ID_HOG (18-char string, unique) as in OG's do-file; the
# destring'd hh_id exceeds double precision and collides (5.07M unique of 5.19M).
ecd <- as.data.table(haven::read_dta(
  file.path(EC_DATA, "CPV2022", "CPV_hh_dwelling.dta"),
  col_select = c("ID_HOG", "aindex", "aindex5")))
ecp <- as.data.table(haven::read_dta(
  file.path(EC_DATA, "CPV2022", "cpv2022_hh_pop.dta"),
  col_select = c("ID_HOG", "total_teen_births", "earliest_teen_birth_yr",
                 "latest_teen_birth_yr")))
num <- function(d) { for (v in setdiff(names(d), "ID_HOG")) set(d, j = v, value = as.numeric(d[[v]])); d }
ecd <- num(ecd); ecp <- num(ecp)
stopifnot(nrow(ecd) == nrow(ecp), uniqueN(ecd$ID_HOG) == nrow(ecd))
ec <- merge(ecd, ecp, by = "ID_HOG")
stopifnot(nrow(ec) == nrow(ecd))
ec[, ID_HOG := NULL]
ec[, any_teen_births := as.integer(total_teen_births > 0)]
rm(ecd, ecp)
ec <- ec[!is.na(aindex5)]
cat(sprintf("EC households: %s\n", format(nrow(ec), big.mark = ",")))

# --- replication check vs OG: reg any_teen_births aindex5 -> -.0562649 -------
m0 <- lm(any_teen_births ~ aindex5, data = ec)
cat(sprintf("EC replication: slope any_teen_births~aindex5 = %+.6f [OG: -0.056265], N = %s [OG: 5,188,827]\n",
            coef(m0)[2], format(nobs(m0), big.mark = ",")))
# --- replicate his post-filter Erreygers (validation of the bridge) ----------
his_post <- ec[latest_teen_birth_yr > 2020 | is.na(latest_teen_birth_yr)]
cat(sprintf("EC replication: Erreygers CI on OG's post filter = %+.5f [OG: -0.02534]\n",
            err_ci(his_post$any_teen_births, his_post$aindex)))

# --- corrected windowed outcomes ---------------------------------------------
ec[, any_post := as.integer(!is.na(latest_teen_birth_yr) & latest_teen_birth_yr >= 2021)]
ec[, any_pre := as.integer(
     (!is.na(earliest_teen_birth_yr) & earliest_teen_birth_yr %between% c(2016, 2020)) |
     (!is.na(latest_teen_birth_yr)   & latest_teen_birth_yr   %between% c(2016, 2020)))]
ec[, any_pre2 := as.integer(
     (!is.na(earliest_teen_birth_yr) & earliest_teen_birth_yr %between% c(2019, 2020)) |
     (!is.na(latest_teen_birth_yr)   & latest_teen_birth_yr   %between% c(2019, 2020)))]
amb <- ec[total_teen_births >= 3 & earliest_teen_birth_yr < 2016 & latest_teen_birth_yr > 2020, .N]
cat(sprintf("EC ambiguous straddlers (>=3 teen births, span across window): %d (%.3f%% of HHs)\n",
            amb, 100 * amb / nrow(ec)))

res_ec <- ec[, .(pre5 = mean(any_pre), pre = mean(any_pre2), post = mean(any_post)), by = aindex5][order(aindex5)]
cat("EC share of HHs with a teen first birth, by quintile:\n")
print(res_ec)
# gradient slopes (per quintile), annualized shares for comparability
sl_pre  <- lm(I(any_pre2 / 2 * 100)  ~ aindex5, data = ec)
sl_post <- lm(I(any_post / 2 * 100) ~ aindex5, data = ec)
cat(sprintf("EC 5y-pre annualized slope (2016-20, reference): %+.4f\n",
            coef(lm(I(any_pre / 5 * 100) ~ aindex5, data = ec))[2]))
long <- rbind(ec[, .(aindex5, y = any_pre2 / 2 * 100, post = 0)],
              ec[, .(aindex5, y = any_post / 2 * 100, post = 1)])
sl_int <- lm(y ~ aindex5 * post, data = long)
cat(sprintf("EC annualized gradient 2y windows (pp/quintile/yr): pre(2019-20) %+.4f, post(2021-22) %+.4f, change %+.4f (SE %.4f)\n",
            coef(sl_pre)[2], coef(sl_post)[2], coef(sl_int)["aindex5:post"],
            sqrt(diag(vcov(sl_int)))["aindex5:post"]))
cat(sprintf("EC Erreygers CI (corrected windows): pre %+.5f | post %+.5f\n",
            err_ci(ec$any_pre2, ec$aindex), err_ci(ec$any_post, ec$aindex)))
ec_plot <- melt(res_ec, id.vars = "aindex5")[, country := "Panel A. Ecuador"]  # EC first (MM 2026-08-03)
rm(ec, his_post, long); gc()

# ============================ DOMINICAN REPUBLIC =============================
POPDIR <- DR_CENSUS
assetcols <- paste0("P15", LETTERS[1:18])
dwellcols <- c("P01","P03","P04","P05","P06","P10","P11","P13","P14","P17","P18","P19")
sel <- c("PROVINCIA","MUNICIPIO","P25_ORDEN","P26_SEXO","P27_EDAD","P65","P67_ANO",
         assetcols, dwellcols)
cat("\nreading DR 2022 census (selected cols)...\n")
d <- fread(file.path(POPDIR, "BD_FINAL_VIVIENDA_HOGAR_PERSONA_XCNPV_PUB.csv"),
           select = sel, showProgress = FALSE)
d[, hid := cumsum(P25_ORDEN == 1L)]
d[, hh_size := .N, by = hid]
cat(sprintf("DR persons: %s | households: %s\n",
            format(nrow(d), big.mark = ","), format(uniqueN(d$hid), big.mark = ",")))

# --- HH wealth index, dissertation v2 recipe (29 items, drop P15A/P15Q) ------
hh <- d[P25_ORDEN == 1L]
for (c in assetcols) hh[, (c) := fifelse(get(c) == 1L, 1, fifelse(get(c) == 2L, 0, NA_real_))]
na9 <- function(x) fifelse(x == 9L, NA_integer_, x)
hh[, `:=`(
  apartment      = as.numeric(P01 == 2L),
  marginal_dwell = as.numeric(P01 %in% c(5L, 6L, 8L)),
  wall_block     = as.numeric(na9(P03) == 1L),
  roof_concrete  = as.numeric(na9(P04) == 1L),
  floor_quality  = as.numeric(P05 %in% c(1L, 3L, 4L, 5L)),
  kitchen_inside = as.numeric(P06 == 1L),
  water_piped_in = as.numeric(P10 == 1L),
  toilet_flush   = as.numeric(na9(P11) == 1L),
  sewer          = fifelse(is.na(na9(P11)), NA_real_,
                           as.numeric(na9(P11) == 1L & !is.na(P13) & P13 == 1L)),
  garbage_coll   = as.numeric(P14 %in% c(1L, 2L)),
  cook_gas_elec  = as.numeric(na9(P18) %in% c(1L, 4L)),
  light_grid     = as.numeric(na9(P19) == 1L),
  crowded        = fifelse(P17 >= 90L, NA_real_,
                           as.numeric(fifelse(P17 == 0L, TRUE, hh_size / P17 > 2)))
)]
pcacols <- c(setdiff(assetcols, c("P15Q", "P15A")),
             c("crowded","apartment","marginal_dwell","wall_block","roof_concrete",
               "floor_quality","kitchen_inside","water_piped_in","toilet_flush","sewer",
               "garbage_coll","cook_gas_elec","light_grid"))
A <- as.matrix(hh[, ..pcacols]); cc <- complete.cases(A)
pc <- prcomp(A[cc, ], center = TRUE, scale. = TRUE)
sc <- rep(NA_real_, nrow(hh)); sc[cc] <- pc$x[, 1]
if (cor(sc[cc], rowSums(A[cc, seq_len(16)])) < 0) sc <- -sc
hh[, wealth := sc]
hh <- hh[!is.na(wealth), .(hid, wealth)]
hh[, q5 := cut(wealth, quantile(wealth, seq(0, 1, .2)), labels = FALSE, include.lowest = TRUE)]
cat(sprintf("DR HHs with wealth index: %s\n", format(nrow(hh), big.mark = ",")))

# --- teen births from LAST live birth (P67) ----------------------------------
w <- d[P26_SEXO == 2L & !is.na(P67_ANO) & P67_ANO %between% c(2010, 2022)]
w[, age_at_birth := P27_EDAD - (2022L - P67_ANO)]
w[, teen_birth := as.integer(age_at_birth %between% c(15, 19))]
tb <- w[teen_birth == 1, .(
  any_pre5 = as.integer(any(P67_ANO %between% c(2016, 2020))),
  any_pre  = as.integer(any(P67_ANO %between% c(2019, 2020))),
  any_post = as.integer(any(P67_ANO %between% c(2021, 2022)))), by = hid]
hh <- merge(hh, tb, by = "hid", all.x = TRUE)
hh[is.na(any_pre),  any_pre  := 0L]
hh[is.na(any_pre5), any_pre5 := 0L]
hh[is.na(any_post), any_post := 0L]

res_dr <- hh[, .(pre5 = mean(any_pre5), pre = mean(any_pre), post = mean(any_post)), by = q5][order(q5)]
cat("DR share of HHs with a teen (last) birth, by quintile:\n")
print(res_dr)
sl_pre_d  <- lm(I(any_pre / 2 * 100)  ~ q5, data = hh)
sl_post_d <- lm(I(any_post / 2 * 100) ~ q5, data = hh)
longd <- rbind(hh[, .(q5, y = any_pre / 2 * 100, post = 0)],
               hh[, .(q5, y = any_post / 2 * 100, post = 1)])
sl_int_d <- lm(y ~ q5 * post, data = longd)
cat(sprintf("DR annualized gradient 2y windows (pp/quintile/yr): pre(2019-20) %+.4f, post(2021-22) %+.4f, change %+.4f (SE %.4f)\n",
            coef(sl_pre_d)[2], coef(sl_post_d)[2], coef(sl_int_d)["q5:post"],
            sqrt(diag(vcov(sl_int_d)))["q5:post"]))
cat(sprintf("DR Erreygers CI: pre %+.5f | post %+.5f\n",
            err_ci(hh$any_pre, hh$wealth), err_ci(hh$any_post, hh$wealth)))
cat("DR caveat: pre-window teen births are undercounted where a later birth\n")
cat("overwrote them (census records only the LAST live birth).\n")
dr_plot <- melt(res_dr, id.vars = "q5")[, country := "Panel B. Dominican Republic"]
setnames(dr_plot, "q5", "aindex5")
sink()

# ============================ FIGURE 5 =======================================
# Curves with visible gradients (MM 2026-07-30): pre = FULL window 2016-2020
# (annualized /5), post = 2021-2022 (annualized /2), both countries.
# DR caveat: last-birth censoring erodes the DR pre window (stated in caption).
pd <- rbind(dr_plot, ec_plot)
pd <- pd[variable %in% c("pre5", "post")]
pd[, pct_yr := fifelse(variable == "pre5", value / 5 * 100, value / 2 * 100)]
pd[, period := factor(fifelse(variable == "pre5", "Pre (2016-2020)", "Post (2021-2022)"),
                      levels = c("Pre (2016-2020)", "Post (2021-2022)"))]
cols <- c("Panel B. Dominican Republic" = COL_POOR, "Panel A. Ecuador" = COL_RICH)
band5 <- dcast(pd, country + aindex5 ~ period, value.var = "pct_yr")
setnames(band5, c("country", "aindex5", "y_pre", "y_post"))
fig5 <- ggplot(pd, aes(aindex5, pct_yr, color = country, linetype = period, alpha = period)) +
  geom_ribbon(data = band5, aes(x = aindex5, ymin = pmin(y_pre, y_post), ymax = pmax(y_pre, y_post), fill = country),
              inherit.aes = FALSE, alpha = .22, show.legend = FALSE) +
  scale_fill_manual(values = cols) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:5, labels = QLAB) +
  scale_color_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("Pre (2016-2020)" = "22", "Post (2021-2022)" = "solid"), name = NULL) +
  scale_alpha_manual(values = c("Pre (2016-2020)" = .5, "Post (2021-2022)" = 1), guide = "none") +
  facet_wrap(~country) +
  guides(linetype = guide_legend(override.aes = list(color = "grey30"))) +
  labs(x = "Household asset-index quintile (2022 census)",
       y = "Households with a teen birth, % per year",
       title = "The household wealth gradient of teen births, before and after COVID-19",
       subtitle = paste0("2022 census microdata, annualized window shares. Ecuador: first births (age at first birth).\n",
                         "DR: last births only, so the pre-window is undercounted where a later birth followed (pre/post change not identified)")) +
  theme_fert()
ggsave(file.path(OUT_FIG, "fig5_census_hh_gradient.pdf"), fig5, width = 9.5, height = 5.2)
ggsave(file.path(OUT_FIG, "fig5_census_hh_gradient.png"), fig5, width = 9.5, height = 5.2, dpi = 200)
saveRDS(list(ec = res_ec, dr = res_dr), file.path(OUT, "census_gradient.rds"))
cat("Saved fig5_census_hh_gradient (curves, full pre window) and log.\n")
