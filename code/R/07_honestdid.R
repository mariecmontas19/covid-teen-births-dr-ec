# 07_honestdid.R — Rambachan & Roth (2023) sensitivity bounds for the DR
# headline DDD event studies (log_dif, with covariates). Relative-magnitudes
# restriction: post-period violations of parallel trends bounded by Mbar times
# the largest pre-period violation. Target parameter: average post-2020 effect.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))
library(HonestDiD)

mods <- readRDS(file.path(OUT, "dr_models_main.rds"))
PRE  <- as.character(2016:2019)
POST <- as.character(2021:2025)

run_honest <- function(m, l_vec) {
  terms <- sprintf("year::%s:poorest", c(PRE, POST))
  beta  <- coef(m)[terms]
  sig   <- vcov(m)[terms, terms]
  stopifnot(!anyNA(beta))
  orig  <- constructOriginalCS(betahat = beta, sigma = sig,
                               numPrePeriods = length(PRE), numPostPeriods = length(POST),
                               l_vec = l_vec)
  sens  <- createSensitivityResults_relativeMagnitudes(
    betahat = beta, sigma = sig,
    numPrePeriods = length(PRE), numPostPeriods = length(POST),
    l_vec = l_vec, Mbarvec = seq(0.5, 2, by = 0.5))
  list(orig = orig, sens = sens)
}
# Two target parameters: (i) average post-2020 effect (2021-2025);
# (ii) average 2023-2025 effect — the window where the DR effect operates
# (2021-22 conceptions were still in the acute-pandemic bust period).
L_ALL  <- rep(1/5, 5)
L_LATE <- c(0, 0, 1/3, 1/3, 1/3)

res <- list()
sink(file.path(OUT, "logs", "07_honestdid.txt"))
for (k in 1:4) {
  for (tgt in c("all", "late")) {
    lv <- if (tgt == "all") L_ALL else L_LATE
    cat(sprintf("\n== log_dif%d: teen vs %s — %s ==\n", k, AG_LABS[AG_COMP[k]],
                if (tgt == "all") "avg 2021-2025 effect" else "avg 2023-2025 effect"))
    r <- run_honest(mods[[paste0("es_logdif", k, "_X")]], lv)
    cat(sprintf("Original 95%% CI: [%.3f, %.3f]\n", r$orig$lb, r$orig$ub))
    print(as.data.frame(r$sens))
    res[[paste0("logdif", k, "_", tgt)]] <- r
  }
}
sink()
# ---- Ecuador (10 pre-periods 2010-2019, 4 post 2021-2024) -------------------
sink(file.path(OUT, "logs", "07_honestdid.txt"), append = TRUE)
ec_mods <- readRDS(file.path(OUT, "ec_models_poorest.rds"))
PRE  <- as.character(2010:2019); POST <- as.character(2021:2024)
L_ALL <- rep(1/4, 4)
for (k in 1:4) {
  cat(sprintf("\n== ECUADOR log_dif%d: teen vs %s — avg 2021-2024 effect ==\n", k, AG_LABS[AG_COMP[k]]))
  r <- run_honest(ec_mods[[paste0("es_logdif", k, "_X")]], L_ALL)
  cat(sprintf("Original 95%% CI: [%.3f, %.3f]\n", r$orig$lb, r$orig$ub))
  print(as.data.frame(r$sens))
  res[[paste0("ec_logdif", k, "_all")]] <- r
}
sink()

# ---- Smoothness-restriction (Delta^SD) bounds --------------------------------
# Appropriate here: pre-trends are FLAT but NOISY. Delta^SD allows post-period
# trend deviations up to M per period (linear-extrapolation violations), rather
# than benchmarking against the largest (noisy) pre-period violation.
run_sd <- function(m, pre, post, l_vec, Mvec = seq(0, .05, .01)) {
  terms <- sprintf("year::%s:poorest", c(pre, post))
  beta <- coef(m)[terms]; sig <- vcov(m)[terms, terms]
  createSensitivityResults(betahat = beta, sigma = sig,
                           numPrePeriods = length(pre), numPostPeriods = length(post),
                           l_vec = l_vec, Mvec = Mvec)
}
sink(file.path(OUT, "logs", "07_honestdid.txt"), append = TRUE)
cat("\n\n================ SMOOTHNESS (Delta^SD) BOUNDS ================\n")
sd_res <- list()
for (k in 1:4) {
  cat(sprintf("\n== DR log_dif%d (vs %s), avg 2023-2025 effect ==\n", k, AG_LABS[AG_COMP[k]]))
  s <- run_sd(mods[[paste0("es_logdif", k, "_X")]], as.character(2016:2019),
              as.character(2021:2025), c(0, 0, 1/3, 1/3, 1/3))
  print(as.data.frame(s)); sd_res[[paste0("dr", k)]] <- s
}
for (k in 1:4) {
  cat(sprintf("\n== EC log_dif%d (vs %s), avg 2021-2024 effect ==\n", k, AG_LABS[AG_COMP[k]]))
  s <- run_sd(ec_mods[[paste0("es_logdif", k, "_X")]], as.character(2010:2019),
              as.character(2021:2024), rep(1/4, 4))
  print(as.data.frame(s)); sd_res[[paste0("ec", k)]] <- s
}
sink()
saveRDS(sd_res, file.path(OUT, "honestdid_smoothness.rds"))
for (nm in names(sd_res)) {
  s <- as.data.frame(sd_res[[nm]])
  ok <- s$M[s$lb > 0]
  cat(sprintf("SD %s: breakdown M = %s (CI excludes 0 up to this per-period trend violation)\n",
              nm, if (length(ok)) sprintf("%.2f", max(ok)) else "none (0 not excluded even at M=0)"))
}
saveRDS(res, file.path(OUT, "dr_honestdid.rds"))

# Console summary: breakdown Mbar (largest Mbar with CI excluding 0)
for (nm in names(res)) {
  s <- as.data.frame(res[[nm]]$sens)
  ok <- s$Mbar[s$lb > 0]
  cat(sprintf("%s: original CI [%.3f, %.3f]; robust to Mbar <= %s\n",
              nm, res[[nm]]$orig$lb, res[[nm]]$orig$ub,
              if (length(ok)) max(ok) else "none"))
}
