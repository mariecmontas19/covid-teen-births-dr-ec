# 05_figures.R — comparative event-study figures (DR + Ecuador), shared theme.
# Fig 1: raw teen birth rates, poorest quintile vs rest, both countries.
# Fig 2: headline log-difference DDD event studies (teen vs each comparison),
#        year x poorest coefficients with 95% CIs, base 2020.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

dr_mods <- readRDS(file.path(OUT, "dr_models_main.rds"))
ec_mods <- readRDS(file.path(OUT, "ec_models_poorest.rds"))
dr_pan  <- as.data.table(readRDS(file.path(OUT, "dr_panel.rds")))
ec_pan  <- as.data.table(haven::read_dta(file.path(EC_DATA, "base_canton_ENV_FE_rev2026.dta"),
                                         col_select = c("canton", "year", "tn1519f", "pob1519f", "poorest")))
ec_pan <- ec_pan[, lapply(.SD, as.numeric)]   # strip haven_labelled attrs

grab <- function(m, country, comp) {
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  ct <- ct[grepl("^year::\\d{4}:poorest$", term)]
  ct[, year := as.integer(gsub("year::(\\d{4}):poorest", "\\1", term))]
  base <- data.table(term = "base", Estimate = 0, `Std. Error` = 0, year = 2020L)
  out <- rbind(ct[, .(term, Estimate, se = `Std. Error`, year)],
               base[, .(term, Estimate, se = `Std. Error`, year)], fill = TRUE)
  out[, `:=`(country = country, comp = AG_LABS[comp],
             lo = Estimate - 1.96 * se, hi = Estimate + 1.96 * se)]
  out
}

es <- rbind(
  rbindlist(lapply(seq_along(AG_COMP), \(k)
    grab(dr_mods[[paste0("es_logdif", k)]], "Dominican Republic", AG_COMP[k]))),
  rbindlist(lapply(seq_along(AG_COMP), \(k)
    grab(ec_mods[[paste0("es_logdif", k)]], "Ecuador", AG_COMP[k]))))
es[, comp := factor(comp, levels = AG_LABS[AG_COMP])]
es[, country := factor(country, levels = c("Ecuador", "Dominican Republic"))]  # EC first (MM 2026-08-03)

# One event-study figure per country (Fig 2 DR, Fig 2 EC)
plot_es <- function(cty, col, unit, units_pl) {
  d <- es[country == cty]
  ggplot(d, aes(year, Estimate)) +
    geom_hline(yintercept = 0, linewidth = .3) +
    geom_vline(xintercept = 2020.5, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(ymin = lo, ymax = hi), color = col, size = .3) +
    scale_x_continuous(breaks = seq(min(d$year), max(d$year), by = 2)) +
    facet_wrap(~comp, ncol = 2,
               labeller = labeller(comp = \(x) paste("Teen (15-19) vs", x))) +
    labs(x = NULL, y = "Year x poorest-quintile coefficient (log difference)",
         title = sprintf("%s: adolescent-to-adult birth-rate differences by local wealth", cty),
         subtitle = sprintf("DDD event studies: log teen rate minus log comparison-group rate; base year 2020; 95%% CIs, SEs clustered by %s", unit)) +
    theme_fert()
}
ggsave(file.path(OUT_FIG, "fig2_DR_eventstudy_logdif_poorest.pdf"),
       plot_es("Dominican Republic", COL_POOR, "municipality", "municipalities"), width = 9, height = 6.5)
ggsave(file.path(OUT_FIG, "fig2_DR_eventstudy_logdif_poorest.png"),
       plot_es("Dominican Republic", COL_POOR, "municipality", "municipalities"), width = 9, height = 6.5, dpi = 200)
ggsave(file.path(OUT_FIG, "fig2_EC_eventstudy_logdif_poorest.pdf"),
       plot_es("Ecuador", COL_RICH, "canton", "cantons"), width = 9, height = 6.5)
ggsave(file.path(OUT_FIG, "fig2_EC_eventstudy_logdif_poorest.png"),
       plot_es("Ecuador", COL_RICH, "canton", "cantons"), width = 9, height = 6.5, dpi = 200)

# Combined figure: both countries in one event-study plot
fig2_both <- ggplot(es, aes(year, Estimate)) +
  geom_hline(yintercept = 0, linewidth = .3) +
  geom_vline(xintercept = 2020.5, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi, color = country),
                  position = position_dodge(width = .45), size = .25) +
  scale_color_manual(values = c("Dominican Republic" = COL_POOR, "Ecuador" = COL_RICH), name = NULL) +
  scale_x_continuous(breaks = seq(2010, 2025, by = 3)) +
  facet_wrap(~comp, ncol = 2,
             labeller = labeller(comp = \(x) paste("Teen (15-19) vs", x))) +
  labs(x = NULL, y = "Year x poorest-quintile coefficient (log difference)",
       title = "Adolescent-to-adult birth-rate differences by local wealth",
       subtitle = "DDD event studies: log teen rate minus log comparison-group rate; base year 2020; 95% CIs, SEs clustered by municipality/canton") +
  theme_fert()
ggsave(file.path(OUT_FIG, "fig2_BOTH_eventstudy_logdif_poorest.pdf"), fig2_both, width = 9.5, height = 7)
ggsave(file.path(OUT_FIG, "fig2_BOTH_eventstudy_logdif_poorest.png"), fig2_both, width = 9.5, height = 7, dpi = 200)

# ---- Fig 1: raw teen rates by group ----------------------------------------
tr <- rbind(
  dr_pan[, .(rate = weighted.mean(tn1519, womenA_15_19, na.rm = TRUE)),
         by = .(year, poorest)][, country := "Dominican Republic"],
  ec_pan[, .(rate = weighted.mean(tn1519f, pob1519f, na.rm = TRUE)),
         by = .(year, poorest)][, country := "Ecuador"])
tr[, grp := fifelse(poorest == 1, "Poorest 20%", "Other 80%")]
tr[, country := factor(country, levels = c("Ecuador", "Dominican Republic"))]  # EC first (MM 2026-08-03)

# One trends figure per country (Fig 1 DR, Fig 1 EC)
plot_tr <- function(cty, tag, census) {
  ggplot(tr[country == cty], aes(year, rate, color = grp)) +
    geom_vline(xintercept = 2020.5, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = .8) + geom_point(size = 1.4) +
    scale_color_manual(values = c("Poorest 20%" = COL_POOR, "Other 80%" = COL_RICH), name = NULL) +
    scale_x_continuous(breaks = scales::breaks_width(2), labels = \(x) sprintf("%d", x)) +
    labs(x = NULL, y = "Teen births per 1,000 women 15-19",
         title = sprintf("%s: teen birth rates by local wealth", cty),
         subtitle = sprintf("Population-weighted means; poorest quintile by %s census asset index", census)) +
    theme_fert()
}
ggsave(file.path(OUT_FIG, "fig1_DR_teen_rates_by_wealth.pdf"),
       plot_tr("Dominican Republic", "DR", "2010"), width = 7, height = 4.6)
ggsave(file.path(OUT_FIG, "fig1_DR_teen_rates_by_wealth.png"),
       plot_tr("Dominican Republic", "DR", "2010"), width = 7, height = 4.6, dpi = 200)
ggsave(file.path(OUT_FIG, "fig1_EC_teen_rates_by_wealth.pdf"),
       plot_tr("Ecuador", "EC", "2010"), width = 7, height = 4.6)
ggsave(file.path(OUT_FIG, "fig1_EC_teen_rates_by_wealth.png"),
       plot_tr("Ecuador", "EC", "2010"), width = 7, height = 4.6, dpi = 200)

fig1_both <- ggplot(tr, aes(year, rate, color = grp)) +
  geom_vline(xintercept = 2020.5, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = .8) + geom_point(size = 1.4) +
  scale_color_manual(values = c("Poorest 20%" = COL_POOR, "Other 80%" = COL_RICH), name = NULL) +
  scale_x_continuous(breaks = scales::breaks_width(2), labels = \(x) sprintf("%d", x)) +
  facet_wrap(~country, scales = "free_x") +
  labs(x = NULL, y = "Teen births per 1,000 women 15-19",
       title = "Teen birth rates by local wealth, before and after COVID-19",
       subtitle = "Population-weighted means; poorest quintile by 2010 census asset index (both countries)") +
  theme_fert()
ggsave(file.path(OUT_FIG, "fig1_BOTH_teen_rates_by_wealth.pdf"), fig1_both, width = 9.5, height = 4.6)
ggsave(file.path(OUT_FIG, "fig1_BOTH_teen_rates_by_wealth.png"), fig1_both, width = 9.5, height = 4.6, dpi = 200)
cat("Figures written to", OUT_FIG, "\n")
