# 02_dr_did_eventstudy.R — DR DiD and event studies mirroring Omar's Ecuador
# specs (xtreg ..., fe vce(cluster canton)), with the log-difference DDD as
# the headline outcome. Base year 2020; post = 2021+.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

panel <- as.data.table(readRDS(file.path(OUT, "dr_panel.rds")))
XVARS <- c("avghaitian", "avgage", "avgeduc", "sr", "avginst", "avgsba", "avgcp")
XRHS  <- paste(XVARS, collapse = " + ")

es  <- function(y, X = FALSE) feols(
  as.formula(sprintf("%s ~ i(year, poorest, ref = 2020) + factor(year)%s | id",
                     y, if (X) paste0(" + ", XRHS) else "")),
  data = panel, cluster = ~id)
did <- function(y, X = FALSE) feols(
  as.formula(sprintf("%s ~ post:poorest + factor(year)%s | id",
                     y, if (X) paste0(" + ", XRHS) else "")),
  data = panel, cluster = ~id)

mods <- list()
# (1) Headline: log-difference DDD event studies (teen vs each comparison group)
for (k in 1:4) {
  mods[[paste0("es_logdif", k)]]   <- es(paste0("log_dif", k))
  mods[[paste0("es_logdif", k, "_X")]] <- es(paste0("log_dif", k), X = TRUE)
}
# (2) Event studies on the log teen rate itself
mods[["es_logteen"]]   <- es("log_tn1519")
mods[["es_logteen_X"]] <- es("log_tn1519", X = TRUE)
# (3) Simple pre/post DiD in logs
mods[["did_logteen"]]    <- did("log_tn1519")
mods[["did_logteen_X"]]  <- did("log_tn1519", X = TRUE)
for (k in 1:4) {
  mods[[paste0("did_logdif", k)]] <- did(paste0("log_dif", k))
  mods[[paste0("did_logdif", k, "_X")]] <- did(paste0("log_dif", k), X = TRUE)
}
# (4) Level-difference event studies (Omar's tn_dif table 2 analog, with X)
for (k in 1:4) mods[[paste0("es_tndif", k, "_X")]] <- es(paste0("tn_dif", k), X = TRUE)

saveRDS(mods, file.path(OUT, "dr_models_main.rds"))

sink(file.path(OUT, "logs", "02_dr_did_eventstudy.txt"))
cat("== DR DiD post x poorest, log teen rate (no X | with X) ==\n")
print(etable(mods$did_logteen, mods$did_logteen_X, se.below = TRUE))
cat("\n== DR DiD post x poorest, log differences (with X) ==\n")
print(etable(mods$did_logdif1_X, mods$did_logdif2_X, mods$did_logdif3_X, mods$did_logdif4_X, se.below = TRUE))
cat("\n== Headline event studies: log_dif (teen - comparison), with X ==\n")
print(etable(mods$es_logdif1_X, mods$es_logdif2_X, mods$es_logdif3_X, mods$es_logdif4_X, se.below = TRUE))
cat("\n== Event study: log teen rate, with X ==\n")
print(etable(mods$es_logteen, mods$es_logteen_X, se.below = TRUE))
sink()
cat("Models saved; summary log written to analysis/output/logs/02_dr_did_eventstudy.txt\n")

# Console teaser: post-2020 coefficients of the headline models
for (k in 1:4) {
  cf <- coef(mods[[paste0("es_logdif", k, "_X")]])
  post_cf <- cf[grep("year::202[1-5]:poorest", names(cf))]
  cat(sprintf("log_dif%d (teen vs %s) year x poorest 2021-2025: %s\n",
              k, AG_LABS[AG_COMP[k]],
              paste(sprintf("%.3f", post_cf), collapse = " ")))
}
