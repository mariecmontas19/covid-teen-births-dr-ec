# 08_continuous_wealth.R — robustness: continuous poverty gradient instead of
# the poorest-quintile dummy, plus the 2022-index quintile for comparison.
# pov_z = standardized NEGATIVE 2010 wealth index (higher = poorer), so a
# positive interaction = teen rates rise more in poorer municipalities.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

panel <- as.data.table(readRDS(file.path(OUT, "dr_panel.rds")))
XRHS  <- paste(c("avghaitian", "avgage", "avgeduc", "sr", "avginst", "avgsba", "avgcp"), collapse = " + ")
u <- unique(panel[, .(id, wealth_index_2010)])
panel[, pov_z := -(wealth_index_2010 - mean(u$wealth_index_2010)) / sd(u$wealth_index_2010)]

mods <- list()
for (k in 1:4) {
  mods[[paste0("did_cont", k)]] <- feols(
    as.formula(sprintf("log_dif%d ~ post:pov_z + factor(year) + %s | id", k, XRHS)),
    data = panel, cluster = ~id)
  mods[[paste0("es_cont", k)]] <- feols(
    as.formula(sprintf("log_dif%d ~ i(year, pov_z, ref = 2020) + factor(year) + %s | id", k, XRHS)),
    data = panel, cluster = ~id)
  mods[[paste0("did_q22_", k)]] <- feols(
    as.formula(sprintf("log_dif%d ~ post:poorest_2022 + factor(year) + %s | id", k, XRHS)),
    data = panel, cluster = ~id)
}

# ---- Ecuador continuous gradient (CPV2010 asset_index, canton level) --------
ec <- as.data.table(haven::read_dta(file.path(EC_DATA, "base_canton_ENV_FE_rev2026.dta"),
        col_select = c("canton", "year", paste0("log_dif", 1:4), "post",
                       "avgind", "avgage", "avgeduc", "sr", "avginst", "avgsba",
                       "avgcp", "avgchild2")))
ec <- ec[, lapply(.SD, as.numeric)]
ai <- as.data.table(haven::read_dta(file.path(EC_DATA, "CPV2010", "CPV2010_aindex_canton.dta")))
ec <- merge(ec, ai[, .(canton = as.numeric(canton), asset_index)], by = "canton")
ue <- unique(ec[, .(canton, asset_index)])
ec[, pov_z := -(asset_index - mean(ue$asset_index)) / sd(ue$asset_index)]
XRHS_EC <- "avgind + avgage + avgeduc + sr + avginst + avgsba + avgcp + avgchild2"
for (k in 1:4) {
  mods[[paste0("ec_did_cont", k)]] <- feols(
    as.formula(sprintf("log_dif%d ~ post:pov_z + factor(year) + %s | canton", k, XRHS_EC)),
    data = ec, cluster = ~canton)
  mods[[paste0("ec_es_cont", k)]] <- feols(
    as.formula(sprintf("log_dif%d ~ i(year, pov_z, ref = 2020) + factor(year) + %s | canton", k, XRHS_EC)),
    data = ec, cluster = ~canton)
}

saveRDS(mods, file.path(OUT, "dr_models_robustness.rds"))

# Comparative continuous-gradient event-study figure (fig4)
grab <- function(m, country, comp) {
  ct <- as.data.table(coeftable(m), keep.rownames = "term")[grepl("^year::\\d{4}:pov_z$", term)]
  ct[, year := as.integer(gsub("year::(\\d{4}):pov_z", "\\1", term))]
  out <- rbind(ct[, .(Estimate, se = `Std. Error`, year)],
               data.table(Estimate = 0, se = 0, year = 2020L))
  out[, `:=`(country = country, comp = AG_LABS[comp],
             lo = Estimate - 1.96 * se, hi = Estimate + 1.96 * se)]
}
esd <- rbind(
  rbindlist(lapply(1:4, \(k) grab(mods[[paste0("es_cont", k)]], "Dominican Republic", AG_COMP[k]))),
  rbindlist(lapply(1:4, \(k) grab(mods[[paste0("ec_es_cont", k)]], "Ecuador", AG_COMP[k]))))
esd[, comp := factor(comp, levels = AG_LABS[AG_COMP])]
esd[, country := factor(country, levels = c("Ecuador", "Dominican Republic"))]  # EC first (MM 2026-08-03)
fig4 <- ggplot(esd, aes(year, Estimate)) +
  geom_hline(yintercept = 0, linewidth = .3) +
  geom_vline(xintercept = 2020.5, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi, color = country),
                  position = position_dodge(width = .45), size = .25) +
  scale_color_manual(values = c("Dominican Republic" = COL_POOR, "Ecuador" = COL_RICH), name = NULL) +
  scale_x_continuous(breaks = seq(2010, 2025, by = 3)) +
  facet_wrap(~comp, ncol = 2, labeller = labeller(comp = \(x) paste("Teen (15-19) vs", x))) +
  labs(x = NULL, y = "Year x poverty-index (1 SD) coefficient (log difference)",
       title = "Continuous poverty gradient of the post-COVID teen birth rebound",
       subtitle = "Standardized 2010-census poverty index (higher = poorer) replaces the quintile dummy.\nReference year 2020; 95% confidence intervals.") +
  theme_fert()
ggsave(file.path(OUT_FIG, "fig3_BOTH_eventstudy_continuous_gradient.pdf"), fig4, width = 9.5, height = 7)
ggsave(file.path(OUT_FIG, "fig3_BOTH_eventstudy_continuous_gradient.png"), fig4, width = 9.5, height = 7, dpi = 200)
sink(file.path(OUT, "logs", "08_continuous_wealth.txt"))
cat("== Continuous poverty gradient: post x pov_z (1 SD poorer) ==\n")
print(etable(mods$did_cont1, mods$did_cont2, mods$did_cont3, mods$did_cont4,
             keep = "%post|%pov", headers = paste("vs", AG_LABS[AG_COMP]), se.below = TRUE))
cat("\n== 2022-index poorest quintile (comparison with main 2010 definition) ==\n")
print(etable(mods$did_q22_1, mods$did_q22_2, mods$did_q22_3, mods$did_q22_4,
             keep = "%post|%poorest", headers = paste("vs", AG_LABS[AG_COMP]), se.below = TRUE))
cat("\n== Event studies, continuous gradient ==\n")
print(etable(mods$es_cont1, mods$es_cont2, mods$es_cont3, mods$es_cont4,
             keep = "%pov", headers = paste("vs", AG_LABS[AG_COMP]), se.below = TRUE))
sink()

for (k in 1:4) {
  r <- coeftable(mods[[paste0("did_cont", k)]])["post:pov_z", ]
  r22 <- coeftable(mods[[paste0("did_q22_", k)]])["post:poorest_2022", ]
  cat(sprintf("vs %s: post x pov_z %+.4f (p=%.3f) | post x poorest2022 %+.4f (p=%.3f)\n",
              AG_LABS[AG_COMP[k]], r[1], r[4], r22[1], r22[4]))
}
