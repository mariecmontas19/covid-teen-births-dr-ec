# 13_concentration_curves.R — formal gradient analysis requested by OG
# (2026-07-29): health-equity concentration curves and concentration indices
# for adolescent births, before and after the pandemic, both countries.
# Unit = locality ranked by the 2010 census asset index (poorest -> richest),
# weighted by teen women; y = cumulative share of teen births. Concentration
# index via the direct grouped formula, cluster-bootstrap CIs (localities),
# validated against rineq::ci (R analog of Stata's conindex framework).
# Periods: DR pre = 2019-2020 (registration-complete), post = 2021-2025;
#          EC pre = 2010-2020, post = 2021-2024.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))
set.seed(20260729)

dr <- as.data.table(readRDS(file.path(OUT, "dr_panel.rds")))
dr_g <- dr[year >= 2019, .(births = sum(n_births_15_19), women = sum(womenA_15_19),
                           wealth = first(wealth_index_2010),
                           period = fifelse(post == 1, "Post", "Pre")),
           by = .(id, post)]
ec_raw <- as.data.table(haven::read_dta(file.path(EC_DATA, "base_canton_ENV_FE_rev2026.dta"),
          col_select = c("canton", "year", "nv1519f", "pob1519f", "post")))
ec_raw <- ec_raw[, lapply(.SD, as.numeric)]
ai <- as.data.table(haven::read_dta(file.path(EC_DATA, "CPV2010", "CPV2010_aindex_canton.dta")))
ec_raw <- merge(ec_raw, ai[, .(canton = as.numeric(canton), asset_index)], by = "canton")
ec_g <- ec_raw[, .(births = sum(nv1519f, na.rm = TRUE), women = sum(pob1519f),
                   wealth = first(asset_index),
                   period = fifelse(post == 1, "Post", "Pre")),
               by = .(canton, post)]
setnames(ec_g, "canton", "id")

# Concentration index for grouped data: C = sum(w_i * y_i * (2R_i - 1)) / mu,
# w = teen-women share, R = women-weighted fractional rank by wealth (poorest first)
conc_index <- function(d) {
  d <- d[order(wealth)]
  w <- d$women / sum(d$women)
  R <- cumsum(w) - w / 2
  y <- d$births / d$women
  mu <- sum(w * y)
  sum(w * y * (2 * R - 1)) / mu
}
curve_pts <- function(d) {
  d <- d[order(wealth)]
  data.table(cw = c(0, cumsum(d$women) / sum(d$women)),
             cb = c(0, cumsum(d$births) / sum(d$births)))
}
boot_ci <- function(d, B = 2000) {
  ids <- unique(d$id)
  reps <- replicate(B, conc_index(d[id %in% sample(ids, replace = TRUE)][
    , .(births = sum(births), women = sum(women), wealth = first(wealth)), by = id]))
  quantile(reps, c(.025, .975))
}
boot_diff <- function(dpre, dpost, B = 2000) {
  ids <- unique(dpre$id)
  reps <- replicate(B, {
    s <- sample(ids, replace = TRUE)
    agg <- function(d) d[id %in% s][, .(births = sum(births), women = sum(women),
                                        wealth = first(wealth)), by = id]
    conc_index(agg(dpost)) - conc_index(agg(dpre))
  })
  c(quantile(reps, c(.025, .975)), p = mean(sign(reps) != sign(median(reps))) * 2)
}

res <- list(); pts <- list()
for (cty in c("DR", "EC")) {
  g <- if (cty == "DR") dr_g else ec_g
  for (per in c("Pre", "Post")) {
    d <- g[period == per]
    ci_b <- boot_ci(d)
    res[[paste(cty, per)]] <- c(C = conc_index(d), lo = ci_b[1], hi = ci_b[2])
    pts[[paste(cty, per)]] <- curve_pts(d)[, `:=`(country = cty, period = per)]
  }
  dd <- boot_diff(g[period == "Pre"], g[period == "Post"])
  res[[paste(cty, "diff")]] <- c(C = conc_index(g[period == "Post"]) - conc_index(g[period == "Pre"]),
                                 lo = dd[1], hi = dd[2])
}

sink(file.path(OUT, "logs", "13_concentration_curves.txt"))
cat("Concentration indices (negative = concentrated among the poor), cluster bootstrap 95% CI\n")
for (nm in names(res)) cat(sprintf("%-8s C = %+.4f [%+.4f, %+.4f]\n",
                                   nm, res[[nm]][1], res[[nm]][2], res[[nm]][3]))
# Validation vs rineq (locality-level rates weighted by women)
if (requireNamespace("rineq", quietly = TRUE)) {
  for (cty in c("DR", "EC")) {
    g <- if (cty == "DR") dr_g else ec_g
    for (per in c("Pre", "Post")) {
      d <- g[period == per][order(wealth)]
      v <- rineq::ci(ineqvar = d$wealth, outcome = d$births / d$women,
                     weights = d$women, type = "CI", method = "direct")
      cat(sprintf("rineq check %s %s: C = %+.4f (ours %+.4f)\n",
                  cty, per, v$concentration_index, res[[paste(cty, per)]][1]))
    }
  }
}
sink()
system(paste("cat", shQuote(file.path(OUT, "logs", "13_concentration_curves.txt"))))
saveRDS(res, file.path(OUT, "concentration_indices.rds"))

# ---- Figure 4: concentration curves — Panel A DR, Panel B EC; pre vs post ----
# (MM 2026-07-29: country panels, pre and post on the same axes, shaded band
#  between the curves so the pre->post movement is visually obvious.)
pd <- rbindlist(pts)
pd[, country := fifelse(country == "DR", "Panel B. Dominican Republic", "Panel A. Ecuador")]  # EC first (MM 2026-08-03)

# Interpolate both curves onto a common grid to shade the area between them
grid <- seq(0, 1, by = .01)
band <- rbindlist(lapply(unique(pd$country), function(cc) {
  pre  <- pd[country == cc & period == "Pre"]
  post <- pd[country == cc & period == "Post"]
  data.table(country = cc, cw = grid,
             y_pre  = approx(pre$cw,  pre$cb,  xout = grid, ties = "ordered")$y,
             y_post = approx(post$cw, post$cb, xout = grid, ties = "ordered")$y)
}))

ann <- data.table(
  country = rep(c("Panel A. Ecuador", "Panel B. Dominican Republic"), each = 2),
  x = .55, y = c(.30, .18, .30, .18),
  period = rep(c("Pre", "Post"), 2),
  txt = c(sprintf("Pre:  C = %.3f", res[["EC Pre"]][1]),
          sprintf("Post: C = %.3f", res[["EC Post"]][1]),
          sprintf("Pre:  C = %.3f", res[["DR Pre"]][1]),
          sprintf("Post: C = %.3f", res[["DR Post"]][1])))

fig4 <- ggplot() +
  geom_ribbon(data = band, aes(x = cw, ymin = y_pre, ymax = y_post, fill = country),
              alpha = .25, show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey75", linewidth = .4) +
  annotate("text", x = .84, y = .78, label = "Line of equality", color = "grey55", size = 3, angle = 38) +
  geom_line(data = pd, aes(cw, cb, color = country, linetype = period,
                           alpha = period), linewidth = .9) +
  geom_text(data = ann, aes(x, y, label = txt, color = country, alpha = period),
            size = 3.4, hjust = 0, show.legend = FALSE) +
  scale_color_manual(values = c("Panel B. Dominican Republic" = COL_POOR,
                                "Panel A. Ecuador" = COL_RICH), guide = "none") +
  scale_fill_manual(values = c("Panel B. Dominican Republic" = COL_POOR,
                               "Panel A. Ecuador" = COL_RICH)) +
  scale_linetype_manual(values = c(Pre = "22", Post = "solid"), name = NULL,
                        labels = c(Pre = "Pre-pandemic", Post = "Post-pandemic (2021+)")) +
  scale_alpha_manual(values = c(Pre = .55, Post = 1), guide = "none") +
  facet_wrap(~country) +
  coord_equal() +
  guides(linetype = guide_legend(override.aes = list(color = "grey30"))) +
  labs(x = "Cumulative share of teen women, poorest to richest localities",
       y = "Cumulative share of teen births",
       title = "Concentration curves of adolescent births, before and after COVID-19",
       subtitle = paste0("Localities ranked by the 2010 census asset index; shaded band = pre-to-post movement.\n",
                         "C = concentration index (negative = concentrated among the poor)")) +
  theme_fert()
ggsave(file.path(OUT_FIG, "fig4_concentration_curves.pdf"), fig4, width = 9.5, height = 5.4)
ggsave(file.path(OUT_FIG, "fig4_concentration_curves.png"), fig4, width = 9.5, height = 5.4, dpi = 200)
cat("Saved fig4_concentration_curves (country panels, pre/post overlay).\n")
