# 00b_wealth2010.R — DR wealth index from the 2010 census (CENSO2010RD-HOGARES),
# mirroring the 2022-census construction in the dissertation's 05_covariates.R
# and Omar's Ecuador CPV2010 index: household-level PCA on durables + dwelling
# services, PC1 oriented so higher = wealthier, standardized at HH level,
# municipality value = mean of household scores. Item-screening rule as in the
# 2022 v2 index: overall KMO >= .90 and all item MSAs >= .85 (drop offenders).
# 2010 codings verified vs Libro de codigo BD_Hogares_Censo2010:
#   H09A-H09O durables (1 Si / 2 No / 0 doble marca / 9 no declarado)
#   H11 dormitorios; H12 sanitario (1 Inodoro); H13 uso exclusivo (1 solo hogar)
#   H14 basura (1/2 recogida); H15 agua (1 acueducto dentro);
#   H16 combustible (1 gas propano / 4 electricidad); H17 alumbrado (1 red publica)
# Note: walls/roof/floor/kitchen are in the (unavailable) 2010 vivienda file ->
# not included; the 2010 index is durables+services based, like Ecuador's.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_config.R"))
suppressPackageStartupMessages(library(psych))

POPDIR <- DR_CENSUS
DISS_CLEAN <- DR_CLEAN

durcols <- c("H09A_ESTUFA","H09B_NEVERA","H09C_LAVADORA","H09D_TELEVISOR","H09E_RADIO",
             "H09F_CISTERNA","H09G_TINACO","H09H_COMPUTADORA","H09I_INTERNET","H09J_INVERSOR",
             "H09K_PLANTA","H09L_TELEFONO_FIJO","H09M_CELULAR","H09N_AIRE","H09NN_MOTOR",
             "H09O_AUTOMOVIL")
sel <- c("PROVINCIA","MUNICIPIO", durcols,
         "H11_DORMITORIOS","H12_SANITARIO","H13_USO_SANITARIO","H14_BASURA",
         "H15_PROCEDENCIA_AGUA","H16_COMBUSTIBLE","H17_ALUMBRADO","H25C_TOTAL")
cat("reading CENSO2010RD-HOGARES (selected cols)...\n")
hh <- fread(file.path(POPDIR, "CENSO2010RD-HOGARES.csv"), select = sel, showProgress = FALSE)
cat("households:", format(nrow(hh), big.mark = ","), "\n")

for (c in durcols) hh[, (c) := fifelse(get(c) == 1L, 1, fifelse(get(c) == 2L, 0, NA_real_))]
hh[, `:=`(
  crowded        = fifelse(is.na(H11_DORMITORIOS) | H11_DORMITORIOS >= 90L, NA_real_,
                           as.numeric(fifelse(H11_DORMITORIOS == 0L, TRUE,
                                              H25C_TOTAL / H11_DORMITORIOS > 2))),
  toilet_flush   = as.numeric(H12_SANITARIO == 1L),
  toilet_sole    = fifelse(H12_SANITARIO == 3L, 0,
                           as.numeric(!is.na(H13_USO_SANITARIO) & H13_USO_SANITARIO == 1L)),
  garbage_coll   = as.numeric(H14_BASURA %in% c(1L, 2L)),
  water_piped_in = as.numeric(H15_PROCEDENCIA_AGUA == 1L),
  cook_gas_elec  = as.numeric(H16_COMBUSTIBLE %in% c(1L, 4L)),
  light_grid     = as.numeric(H17_ALUMBRADO == 1L)
)]
svc <- c("crowded","toilet_flush","toilet_sole","garbage_coll","water_piped_in",
         "cook_gas_elec","light_grid")

# ---- item screening: KMO >= .90 overall, item MSAs >= .85 (2022 v2 rule) ----
pcacols <- c(durcols, svc)
repeat {
  A  <- as.matrix(hh[, ..pcacols]); cc <- complete.cases(A)
  km <- psych::KMO(cor(A[cc, ]))
  worst <- names(which.min(km$MSAi))
  cat(sprintf("KMO %.4f | min item MSA %.3f (%s) | items %d | complete HH %s\n",
              km$MSA, min(km$MSAi), worst, length(pcacols), format(sum(cc), big.mark = ",")))
  if (min(km$MSAi) >= .85) break
  cat("  -> dropping", worst, "\n")
  pcacols <- setdiff(pcacols, worst)
}
stopifnot(km$MSA >= .90)
# Save the final KMO detail for the appendix KMO table (tabA8; MM 2026-08-03)
saveRDS(list(items = pcacols, msai = km$MSAi, overall = km$MSA,
             n_complete = sum(cc)), file.path(OUT, "wealth2010_kmo.rds"))

pc <- prcomp(A[cc, ], center = TRUE, scale. = TRUE)
sc <- rep(NA_real_, nrow(hh)); sc[cc] <- pc$x[, 1]
if (cor(sc[cc], rowSums(A[cc, intersect(durcols, pcacols)])) < 0) sc <- -sc
hh[, wealth := as.numeric(scale(sc))]
cat(sprintf("PC1 explains %.1f%% | final items: %s\n",
            100 * summary(pc)$importance[2, 1], paste(pcacols, collapse = ", ")))

# ---- aggregate to municipality (2010 codes -> adm3_pcode, 155 munis) --------
m <- as.data.table(readRDS(file.path(DISS_CLEAN, "census_code_to_adm3.rds")))
hh <- merge(hh, m, by.x = c("PROVINCIA", "MUNICIPIO"), by.y = c("prov_code", "muni_code"), all.x = TRUE)
cat("unmapped households:", hh[is.na(adm3_pcode), .N], "\n")
mu <- hh[!is.na(adm3_pcode), .(wealth_index_2010 = mean(wealth, na.rm = TRUE), n_hh_2010 = .N),
         by = adm3_pcode]
# 3 municipios created in 2013 inherit their 2010 parent's value (same territory)
fold <- c(DOM010905 = "DOM010901", DOM051703 = "DOM051701", DOM012510 = "DOM012501")
add <- data.table(adm3_pcode = names(fold),
                  wealth_index_2010 = mu$wealth_index_2010[match(fold, mu$adm3_pcode)],
                  n_hh_2010 = NA_integer_)
mu <- rbind(mu[!adm3_pcode %in% names(fold)], add)

# ---- validation vs 2022 index and SIUBEN ICV poverty ------------------------
cc22 <- as.data.table(readRDS(file.path(DR_CLEAN, "census_covariates_muni.rds")))[, .(adm3_pcode, wealth_index)]
icv  <- as.data.table(readRDS(file.path(DR_CLEAN, "icv_poverty_muni.rds")))[, .(adm3_pcode, pct_poor_icv)]
v <- Reduce(function(a, b) merge(a, b, by = "adm3_pcode"), list(mu, cc22, icv))
cat(sprintf("corr(2010 index, 2022 index) = %.3f | corr(2010 index, ICV poverty) = %.3f | munis matched: %d\n",
            cor(v$wealth_index_2010, v$wealth_index), cor(v$wealth_index_2010, v$pct_poor_icv), nrow(v)))
q20_10 <- quantile(mu$wealth_index_2010, .20)
v[, poorest10 := wealth_index_2010 <= q20_10]
v[, poorest22 := wealth_index <= quantile(wealth_index, .20)]
cat(sprintf("poorest-quintile overlap 2010 vs 2022: %d / %d munis\n",
            v[poorest10 & poorest22, .N], v[poorest10 == TRUE, .N]))

saveRDS(mu, file.path(OUT, "wealth2010_muni.rds"))
cat("Saved:", file.path(OUT, "wealth2010_muni.rds"), "\n")
