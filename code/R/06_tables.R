# 06_tables.R — publication-format .tex tables for the preliminary results.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))

dr  <- readRDS(file.path(OUT, "dr_models_main.rds"))
drp <- readRDS(file.path(OUT, "dr_models_preferred.rds"))
ec  <- readRDS(file.path(OUT, "ec_models_poorest.rds"))

dict <- c(post = "Post (2021+)", poorest = "Poorest quintile",
          "post:poorest" = "Post x Poorest", id = "Municipality", canton = "Canton",
          year = "Year", avghaitian = "Haitian mother share", avgind = "Indigenous mother share",
          avgage = "Avg. mother age", avgeduc = "Avg. mother education", sr = "Sex ratio",
          avginst = "Institutional delivery", avgsba = "Skilled attendant",
          avgcp = "4+ prenatal controls", avgchild2 = "2nd+ order births")
setFixest_dict(dict)
st <- list("-^Fixed-effects", "-^Slopes")

# Table 1: DR simple DiD (log teen rate + the four log differences)
etable(dr$did_logteen_X, dr$did_logdif1_X, dr$did_logdif2_X, dr$did_logdif3_X, dr$did_logdif4_X,
       keep = "%post:poorest", headers = c("log teen", "vs 20-24", "vs 25-29", "vs 30-34", "vs 35-49"),
       tex = TRUE, file = file.path(OUT_TAB, "tab1_dr_did.tex"), replace = TRUE,
       title = "DR: Post-COVID x poorest-quintile DiD (municipality FE, clustered SEs)")

# Table 2: DR headline event studies (log differences, with covariates)
etable(dr$es_logdif1_X, dr$es_logdif2_X, dr$es_logdif3_X, dr$es_logdif4_X,
       keep = "%poorest", headers = c("vs 20-24", "vs 25-29", "vs 30-34", "vs 35-49"),
       tex = TRUE, file = file.path(OUT_TAB, "tab2_dr_eventstudy_logdif.tex"), replace = TRUE,
       title = "DR: DDD event studies, log teen minus log comparison rate (base 2020)")

# Table 3: DR Omar-preferred spec (comparison rate as control)
etable(drp$pref_2024, drp$pref_2529, drp$pref_3034, drp$pref_3549,
       keep = "%poorest|%log_tn", headers = c("ctrl 20-24", "ctrl 25-29", "ctrl 30-34", "ctrl 35-49"),
       tex = TRUE, file = file.path(OUT_TAB, "tab3_dr_pref_ratecontrol.tex"), replace = TRUE,
       title = "DR: Omar's preferred specification (comparison log rate as control)")

# Table 4: Ecuador replication + the missing 30-34 model
etable(ec$pref_2024, ec$pref_2529, ec$pref_3034, ec$pref_3549,
       keep = "%poorest|%log_tn", headers = c("ctrl 20-24", "ctrl 25-29", "ctrl 30-34 (new)", "ctrl 35-49"),
       tex = TRUE, file = file.path(OUT_TAB, "tab4_ec_pref_replicated.tex"), replace = TRUE,
       title = "Ecuador: replication of OG's final regressions (fixest = xtreg)")

cat("Tables written to", OUT_TAB, "\n")
