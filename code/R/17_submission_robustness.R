# ============================================================
# Submission robustness checks: pooled DDD, event studies,
# pretrend tests, adult placebos, zero rates, denominator series,
# and wild-cluster-bootstrap availability
# Authors: Marie Montas and Omar Galarraga
# Purpose: Resolve empirical checks requested during manuscript review.
# Inputs: analysis/output/dr_panel.rds; Ecuador final canton panel;
#         repository-provided Dominican population series A and B.
# Outputs: RDS objects under analysis/output/submission_robustness/;
#          analysis/output/results_summary.md;
#          analysis/output/logs/17_submission_robustness.txt.
# ============================================================

# 0. Setup ----
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript so its location can be resolved reproducibly.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
project_library <- file.path(script_dir, "library")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))
source(file.path(script_dir, "00_config.R"))

ROB_OUT <- file.path(OUT, "submission_robustness")
LOG_FILE <- file.path(OUT, "logs", "17_submission_robustness.txt")
SUMMARY_FILE <- file.path(OUT, "results_summary.md")
dir.create(ROB_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)
on.exit({
  while (sink.number(type = "message") > 2L) sink(type = "message")
  while (sink.number() > 0L) sink()
}, add = TRUE)

cat("17_submission_robustness.R\n")
cat("Started:", format(Sys.time(), tz = "America/Chicago"), "\n")
cat("R:", R.version.string, "\n")
cat("fixest:", as.character(utils::packageVersion("fixest")), "\n\n")

set.seed(20260916)

XR_DR <- c("avghaitian", "avgage", "avgeduc", "sr", "avginst", "avgsba", "avgcp")
XR_EC <- c("avgind", "avgage", "avgeduc", "sr", "avginst", "avgsba", "avgcp", "avgchild2")
SPEC_LABELS <- c(no_covariates = "No covariates", covariates = "Registry covariates")
COMP_LABELS <- unname(AG_LABS[AG_COMP])

stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.10, "*", ""))))
}

fmt_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3L)))
}

write_md_table <- function(x, columns = names(x), digits = 3L) {
  d <- as.data.frame(x[, ..columns])
  for (j in seq_along(d)) {
    if (is.numeric(d[[j]])) {
      if (names(d)[j] %in% c("p", "Asymptotic p", "WCB p")) {
        d[[j]] <- fmt_p(d[[j]])
      } else if (names(d)[j] %in% c(
        "N", "Missing", "Locality-years", "Nonmissing", "Zeros", "localities",
        "years", "first_year", "last_year", "locality_years", "df1", "df2", "Draws"
      )) {
        d[[j]] <- fmt_num(d[[j]], 0L)
      } else {
        d[[j]] <- fmt_num(d[[j]], digits)
      }
    }
    d[[j]][is.na(d[[j]])] <- "NA"
    d[[j]] <- gsub("\\|", "\\\\|", as.character(d[[j]]))
  }
  c(
    paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
    apply(d, 1L, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  )
}

extract_one <- function(model, term) {
  ct <- fixest::coeftable(model)
  if (!term %in% rownames(ct)) {
    return(data.table(
      term = term, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_,
      p_value = NA_real_, n_obs = stats::nobs(model)
    ))
  }
  z <- ct[term, ]
  data.table(
    term = term,
    estimate = unname(z[1L]),
    std_error = unname(z[2L]),
    statistic = unname(z[3L]),
    p_value = unname(z[4L]),
    n_obs = stats::nobs(model)
  )
}

extract_event <- function(model) {
  terms <- grep("^year::[0-9]{4}:poorest$", names(stats::coef(model)), value = TRUE)
  ans <- rbindlist(lapply(terms, function(term) extract_one(model, term)))
  ans[, year := as.integer(sub("^year::([0-9]{4}):poorest$", "\\1", term))]
  setcolorder(ans, c("year", "term", "estimate", "std_error", "statistic", "p_value", "n_obs"))
  ans[]
}

fit_main_models <- function(data, country, id_var, xvars) {
  out <- list(pooled = list(), event = list())
  for (k in seq_along(AG_COMP)) {
    outcome <- paste0("log_dif", k)
    for (spec in names(SPEC_LABELS)) {
      rhs_x <- if (spec == "covariates") paste0(" + ", paste(xvars, collapse = " + ")) else ""
      pooled_fml <- stats::as.formula(sprintf(
        "%s ~ post:poorest + factor(year)%s | %s", outcome, rhs_x, id_var
      ))
      event_fml <- stats::as.formula(sprintf(
        "%s ~ i(year, poorest, ref = 2020) + factor(year)%s | %s",
        outcome, rhs_x, id_var
      ))
      key <- paste(country, paste0("vs_", AG_COMP[k]), spec, sep = "__")
      cluster_fml <- stats::as.formula(paste0("~", id_var))
      out$pooled[[key]] <- fixest::feols(pooled_fml, data = data, cluster = cluster_fml)
      out$event[[key]] <- fixest::feols(event_fml, data = data, cluster = cluster_fml)
    }
  }
  out
}

summarize_pooled <- function(models) {
  rbindlist(lapply(names(models), function(key) {
    parts <- strsplit(key, "__", fixed = TRUE)[[1L]]
    z <- extract_one(models[[key]], "post:poorest")
    z[, `:=`(
      country = parts[1L], comparison = sub("^vs_", "", parts[2L]), specification = parts[3L],
      exact_percent = 100 * (exp(estimate) - 1)
    )]
    z
  }), use.names = TRUE, fill = TRUE)
}

summarize_events <- function(models) {
  rbindlist(lapply(names(models), function(key) {
    parts <- strsplit(key, "__", fixed = TRUE)[[1L]]
    z <- extract_event(models[[key]])
    z[, `:=`(
      country = parts[1L], comparison = sub("^vs_", "", parts[2L]), specification = parts[3L],
      exact_percent = 100 * (exp(estimate) - 1)
    )]
    z
  }), use.names = TRUE, fill = TRUE)
}

joint_pretrend_tests <- function(models, data_by_country) {
  rbindlist(lapply(names(models), function(key) {
    parts <- strsplit(key, "__", fixed = TRUE)[[1L]]
    country <- parts[1L]
    pre_years <- sort(unique(data_by_country[[country]]$year[data_by_country[[country]]$year < 2020]))
    pattern <- sprintf("^year::(%s):poorest$", paste(pre_years, collapse = "|"))
    test <- fixest::wald(models[[key]], keep = pattern, print = FALSE)
    data.table(
      country = country,
      comparison = sub("^vs_", "", parts[2L]),
      specification = parts[3L],
      pre_years = paste(pre_years, collapse = ", "),
      statistic = unname(test$stat),
      p_value = unname(test$p),
      df1 = unname(test$df1),
      df2 = unname(test$df2),
      n_obs = stats::nobs(models[[key]])
    )
  }), use.names = TRUE, fill = TRUE)
}

# 1. Data loading and verification ----
dr_path <- file.path(OUT, "dr_panel.rds")
ec_path <- file.path(EC_DATA, "base_canton_ENV_FE_rev2026.dta")
if (!file.exists(dr_path)) stop("Missing Dominican panel: ", dr_path)
if (!file.exists(ec_path)) stop("Missing Ecuador panel: ", ec_path)

dr <- as.data.table(readRDS(dr_path))
ec_vars <- c(
  "canton", "year", "post", "poorest",
  paste0("tn", c(AG_TEEN, AG_COMP), "f"),
  paste0("log_tn", c(AG_TEEN, AG_COMP), "f"),
  paste0("log_dif", seq_along(AG_COMP)), XR_EC
)
ec <- as.data.table(haven::read_dta(ec_path, col_select = tidyselect::all_of(ec_vars)))
ec <- ec[, lapply(.SD, as.numeric)]

required_dr <- c("id", "year", "post", "poorest", paste0("tn", c(AG_TEEN, AG_COMP)),
                 paste0("log_tn", c(AG_TEEN, AG_COMP)), paste0("log_dif", 1:4), XR_DR)
required_ec <- ec_vars
if (length(setdiff(required_dr, names(dr))) > 0L) {
  stop("Dominican panel missing: ", paste(setdiff(required_dr, names(dr)), collapse = ", "))
}
if (length(setdiff(required_ec, names(ec))) > 0L) {
  stop("Ecuador panel missing: ", paste(setdiff(required_ec, names(ec)), collapse = ", "))
}
stopifnot(uniqueN(dr$id) == 155L, nrow(dr) == 1550L)
stopifnot(uniqueN(ec$canton) == 221L, nrow(ec) == 3315L)

data_inventory <- data.table(
  country = c("DR", "Ecuador"),
  locality_variable = c("id", "canton"),
  localities = c(uniqueN(dr$id), uniqueN(ec$canton)),
  years = c(uniqueN(dr$year), uniqueN(ec$year)),
  first_year = c(min(dr$year), min(ec$year)),
  last_year = c(max(dr$year), max(ec$year)),
  locality_years = c(nrow(dr), nrow(ec))
)
saveRDS(data_inventory, file.path(ROB_OUT, "data_inventory.rds"))
cat("Verified analysis panels:\n")
print(data_inventory)

# The saved comparative panel uses Series A for ages 15--34, but its 35--49
# outcome was assembled separately from the 158-municipality census-anchored
# population file. Rebuild that age group on the same 155-municipality Series A
# geography before fitting the Dominican models.
fold_map <- c(
  DOM010905 = "DOM010901", # San Victor -> Moca
  DOM051703 = "DOM051701", # Baitoa -> Santiago
  DOM012510 = "DOM012501"  # Matanzas -> Bani
)
to_155 <- function(x) data.table::fifelse(x %in% names(fold_map), unname(fold_map[x]), x)

dr_births_path <- file.path(DR_CLEAN, "births_clean_2016_2025.rds")
pop_a_path <- file.path(DR_CLEAN, "pop_municipio_2016_2025_A.rds")
if (!file.exists(dr_births_path)) stop("Missing Dominican births file: ", dr_births_path)
if (!file.exists(pop_a_path)) stop("Missing Dominican Series A population file: ", pop_a_path)

dr_births <- as.data.table(readRDS(dr_births_path))
dr_births[, id_155 := to_155(adm3_pcode)]
b3549_a <- dr_births[
  !is.na(id_155) & !is.na(birth_year) & age_mom >= 35 & age_mom <= 49,
  .(nb_3549_a = .N), by = .(id = id_155, year = birth_year)
]
pop_a <- as.data.table(readRDS(pop_a_path))
w3549_a <- pop_a[
  sex == "F" & age_group %in% c("35-39", "40-44", "45-49"),
  .(women3549_a = sum(pop)), by = .(id = adm3_pcode, year)
]
dr <- merge(dr, b3549_a, by = c("id", "year"), all.x = TRUE)
dr <- merge(dr, w3549_a, by = c("id", "year"), all.x = TRUE)
dr[is.na(nb_3549_a), nb_3549_a := 0L]
stopifnot(nrow(dr) == 1550L, !anyNA(dr$women3549_a), all(dr$women3549_a > 0))
dr_3549_reconstruction_audit <- dr[, .(
  locality_years = .N,
  birth_cells_changed = sum(nb_3549 != nb_3549_a),
  births_original = sum(nb_3549),
  births_corrected = sum(nb_3549_a),
  denominator_cells_changed = sum(abs(women3549 - women3549_a) > 1e-8),
  mean_denominator_percent_change = mean(100 * (women3549_a / women3549 - 1)),
  median_denominator_percent_change = stats::median(100 * (women3549_a / women3549 - 1))
)]
saveRDS(dr_3549_reconstruction_audit, file.path(ROB_OUT, "dr_3549_reconstruction_audit.rds"))
cat("\nDominican 35-49 Series A reconstruction audit:\n")
print(dr_3549_reconstruction_audit)
dr[, `:=`(
  nb_3549 = nb_3549_a,
  women3549 = women3549_a,
  tn3549 = 1000 * nb_3549_a / women3549_a
)]
dr[, log_tn3549 := log(tn3549 + 1)]
dr[, log_dif4 := log_tn1519 - log_tn3549]
dr[, c("nb_3549_a", "women3549_a") := NULL]

# Adult log differences use exactly the headline log(rate + 1) functional form.
dr[, `:=`(
  adult_logdif_2529 = log_tn2024 - log_tn2529,
  adult_logdif_3034 = log_tn2024 - log_tn3034,
  adult_logdif_3549 = log_tn2024 - log_tn3549
)]
ec[, `:=`(
  adult_logdif_2529 = log_tn2024f - log_tn2529f,
  adult_logdif_3034 = log_tn2024f - log_tn3034f,
  adult_logdif_3549 = log_tn2024f - log_tn3549f
)]

# 2. Main pooled DDD and event-study models, with and without covariates ----
dr_models <- fit_main_models(dr, "DR", "id", XR_DR)
ec_models <- fit_main_models(ec, "Ecuador", "canton", XR_EC)
main_models <- list(
  pooled = c(dr_models$pooled, ec_models$pooled),
  event = c(dr_models$event, ec_models$event)
)
saveRDS(main_models, file.path(ROB_OUT, "main_models.rds"))

pooled_results <- summarize_pooled(main_models$pooled)
event_results <- summarize_events(main_models$event)
setorder(pooled_results, country, comparison, specification)
setorder(event_results, country, comparison, specification, year)
saveRDS(pooled_results, file.path(ROB_OUT, "pooled_results.rds"))
saveRDS(event_results, file.path(ROB_OUT, "event_study_results.rds"))

cat("\nMain pooled Post x Poorest estimates:\n")
print(pooled_results[, .(country, comparison, specification, estimate, std_error, p_value, n_obs)])

# 3. Joint tests of every pre-2020 year x poorest coefficient ----
pretrend_tests <- joint_pretrend_tests(
  main_models$event,
  data_by_country = list(DR = dr, Ecuador = ec)
)
setorder(pretrend_tests, country, comparison, specification)
saveRDS(pretrend_tests, file.path(ROB_OUT, "joint_pretrend_tests.rds"))
cat("\nJoint pretrend tests:\n")
print(pretrend_tests)

# 4. Adult log-difference placebo models ----
fit_adult_placebos <- function(data, country, id_var, xvars) {
  models <- list()
  results <- list()
  comps <- c("2529", "3034", "3549")
  for (comp in comps) {
    outcome <- paste0("adult_logdif_", comp)
    for (spec in names(SPEC_LABELS)) {
      rhs_x <- if (spec == "covariates") paste0(" + ", paste(xvars, collapse = " + ")) else ""
      fml <- stats::as.formula(sprintf(
        "%s ~ post:poorest + factor(year)%s | %s", outcome, rhs_x, id_var
      ))
      key <- paste(country, paste0("2024_vs_", comp), spec, sep = "__")
      cluster_fml <- stats::as.formula(paste0("~", id_var))
      models[[key]] <- fixest::feols(fml, data = data, cluster = cluster_fml)
      z <- extract_one(models[[key]], "post:poorest")
      z[, `:=`(
        country = country, comparison = paste0("20-24 vs ", AG_LABS[comp]),
        specification = spec, exact_percent = 100 * (exp(estimate) - 1)
      )]
      results[[key]] <- z
    }
  }
  list(models = models, results = rbindlist(results, use.names = TRUE, fill = TRUE))
}

dr_placebo <- fit_adult_placebos(dr, "DR", "id", XR_DR)
ec_placebo <- fit_adult_placebos(ec, "Ecuador", "canton", XR_EC)
adult_placebo_models <- c(dr_placebo$models, ec_placebo$models)
adult_placebo_results <- rbindlist(list(dr_placebo$results, ec_placebo$results))
setorder(adult_placebo_results, country, comparison, specification)
saveRDS(adult_placebo_models, file.path(ROB_OUT, "adult_placebo_models.rds"))
saveRDS(adult_placebo_results, file.path(ROB_OUT, "adult_placebo_results.rds"))
cat("\nAdult log-difference placebo estimates:\n")
print(adult_placebo_results[, .(country, comparison, specification, estimate, std_error, p_value, n_obs)])

# 5. Summary statistics and prevalence of zero locality-year rates ----
rate_summary_one <- function(data, country, vars, labels) {
  rbindlist(lapply(seq_along(vars), function(j) {
    x <- data[[vars[j]]]
    data.table(
      country = country, age_group = labels[j], variable = vars[j],
      n_nonmissing = sum(!is.na(x)), n_missing = sum(is.na(x)),
      mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE),
      std_dev = stats::sd(x, na.rm = TRUE), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE)
    )
  }))
}

zero_summary_one <- function(data, country, vars, labels) {
  rbindlist(lapply(seq_along(vars), function(j) {
    x <- data[[vars[j]]]
    data.table(
      country = country, age_group = labels[j], variable = vars[j],
      locality_years = length(x), nonmissing = sum(!is.na(x)), missing = sum(is.na(x)),
      zero_rates = sum(x == 0, na.rm = TRUE),
      zero_percent_nonmissing = 100 * mean(x == 0, na.rm = TRUE)
    )
  }))
}

zero_by_wealth_one <- function(data, country, vars, labels) {
  rbindlist(lapply(seq_along(vars), function(j) {
    data[, .(
      locality_years = .N,
      nonmissing = sum(!is.na(get(vars[j]))),
      missing = sum(is.na(get(vars[j]))),
      zero_rates = sum(get(vars[j]) == 0, na.rm = TRUE),
      zero_percent_nonmissing = 100 * mean(get(vars[j]) == 0, na.rm = TRUE)
    ), by = poorest][, `:=`(country = country, age_group = labels[j], variable = vars[j])]
  }))
}

age_labels <- unname(AG_LABS[c(AG_TEEN, AG_COMP)])
dr_rate_vars <- paste0("tn", c(AG_TEEN, AG_COMP))
ec_rate_vars <- paste0("tn", c(AG_TEEN, AG_COMP), "f")
rate_summary <- rbindlist(list(
  rate_summary_one(dr, "DR", dr_rate_vars, age_labels),
  rate_summary_one(ec, "Ecuador", ec_rate_vars, age_labels)
))
zero_rate_summary <- rbindlist(list(
  zero_summary_one(dr, "DR", dr_rate_vars, age_labels),
  zero_summary_one(ec, "Ecuador", ec_rate_vars, age_labels)
))
zero_rate_by_wealth <- rbindlist(list(
  zero_by_wealth_one(dr, "DR", dr_rate_vars, age_labels),
  zero_by_wealth_one(ec, "Ecuador", ec_rate_vars, age_labels)
))
saveRDS(rate_summary, file.path(ROB_OUT, "rate_summary_statistics.rds"))
saveRDS(zero_rate_summary, file.path(ROB_OUT, "zero_rate_summary.rds"))
saveRDS(zero_rate_by_wealth, file.path(ROB_OUT, "zero_rate_by_wealth.rds"))
cat("\nZero-rate prevalence:\n")
print(zero_rate_summary)

# 6. Wild-cluster-bootstrap inference ----
# Use the verified fwildclusterboot implementation for null-imposed WCR11 tests
# of the pooled Post x Poorest coefficient, clustered at the locality level.
wild_candidates <- c("fwildclusterboot", "wildrwolf", "boottest", "WildBootTests")
wild_bootstrap_audit <- data.table(
  package = wild_candidates,
  installed = vapply(wild_candidates, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(wild_candidates, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }, character(1)),
  verified_api = FALSE,
  api = NA_character_,
  status = "Not run: package is not installed"
)

# Verify the API rather than guessing if a recognized package is available.
if (wild_bootstrap_audit[package == "fwildclusterboot", installed]) {
  exports <- getNamespaceExports("fwildclusterboot")
  has_api <- "boottest" %in% exports
  wild_bootstrap_audit[package == "fwildclusterboot", `:=`(
    verified_api = has_api,
    api = if (has_api) "fwildclusterboot::boottest" else NA_character_,
    status = if (has_api) "Available but not automatically run pending model-class verification" else "Installed; expected API not exported"
  )]
}
if (wild_bootstrap_audit[package == "wildrwolf", installed]) {
  exports <- getNamespaceExports("wildrwolf")
  has_api <- "rwolf" %in% exports
  wild_bootstrap_audit[package == "wildrwolf", `:=`(
    verified_api = has_api,
    api = if (has_api) "wildrwolf::rwolf" else NA_character_,
    status = if (has_api) "Available but not automatically run pending model-class verification" else "Installed; expected API not exported"
  )]
}

if (!wild_bootstrap_audit[package == "fwildclusterboot", verified_api]) {
  stop("fwildclusterboot::boottest is required but was not found in the project library.")
}

WILD_B <- 9999L
WILD_SEED <- 20260918L
# fwildclusterboot 0.14.3 requires fixed-effect identifiers to be factors when
# it reconstructs fixest model frames. This changes no fitted values or
# coefficients; it only makes the categorical role explicit for preprocessing.
dr_boot <- copy(dr)
ec_boot <- copy(ec)
dr_boot[, id := factor(id)]
ec_boot[, canton := factor(canton)]
wild_bootstrap_objects <- vector("list", length(main_models$pooled))
names(wild_bootstrap_objects) <- names(main_models$pooled)
wild_bootstrap_rows <- vector("list", length(main_models$pooled))

for (i in seq_along(main_models$pooled)) {
  key <- names(main_models$pooled)[i]
  parts <- strsplit(key, "__", fixed = TRUE)[[1L]]
  cluster_var <- if (parts[1L] == "DR") "id" else "canton"
  comparison_index <- match(sub("^vs_", "", parts[2L]), AG_COMP)
  xvars <- if (parts[1L] == "DR") XR_DR else XR_EC
  rhs_x <- if (parts[3L] == "covariates") paste0(" + ", paste(xvars, collapse = " + ")) else ""
  bootstrap_fml <- stats::as.formula(sprintf(
    "log_dif%s ~ post:poorest + factor(year)%s | %s",
    comparison_index, rhs_x, cluster_var
  ))
  # Fit in the script environment with explicit panel and cluster references.
  # boottest updates the original fixest call internally, so models fitted in
  # the helper function (whose call records `data = data`) are not reusable here.
  bootstrap_model <- if (parts[1L] == "DR") {
    fixest::feols(bootstrap_fml, data = dr_boot)
  } else {
    fixest::feols(bootstrap_fml, data = ec_boot)
  }
  bootstrap_model$call$fml <- bootstrap_fml
  bootstrap_model$call$data <- if (parts[1L] == "DR") quote(dr_boot) else quote(ec_boot)
  bootstrap_model$call_env <- globalenv()
  set.seed(WILD_SEED + i)
  dqrng::dqset.seed(WILD_SEED + i)
  boot_fit <- fwildclusterboot::boottest(
    bootstrap_model,
    param = "post:poorest",
    clustid = cluster_var,
    B = WILD_B,
    type = "rademacher",
    impose_null = TRUE,
    bootstrap_type = "fnw11",
    conf_int = FALSE,
    nthreads = 1L
  )
  wild_bootstrap_objects[[i]] <- boot_fit
  bt <- as.data.table(fwildclusterboot::tidy(boot_fit))
  wild_bootstrap_rows[[i]] <- data.table(
    country = parts[1L],
    comparison = sub("^vs_", "", parts[2L]),
    specification = parts[3L],
    estimate = bt$estimate[1L],
    bootstrap_statistic = bt$statistic[1L],
    bootstrap_p_value = bt$p.value[1L],
    bootstrap_draws = WILD_B,
    bootstrap_type = "WCR11, null imposed, Rademacher weights",
    cluster_variable = cluster_var,
    seed = WILD_SEED + i
  )
}

wild_bootstrap_results <- rbindlist(wild_bootstrap_rows)
wild_bootstrap_results <- merge(
  wild_bootstrap_results,
  pooled_results[, .(country, comparison, specification, asymptotic_p_value = p_value)],
  by = c("country", "comparison", "specification"),
  all.x = TRUE,
  sort = FALSE
)
wild_bootstrap_results[, `:=`(
  asymptotic_significant_5pct = asymptotic_p_value < 0.05,
  bootstrap_significant_5pct = bootstrap_p_value < 0.05,
  significance_changed_5pct = (asymptotic_p_value < 0.05) != (bootstrap_p_value < 0.05)
)]
setorder(wild_bootstrap_results, country, comparison, specification)
wild_bootstrap_audit[package == "fwildclusterboot", status :=
  paste0("Verified and run: ", WILD_B, " null-imposed WCR11 Rademacher draws per pooled model")]
saveRDS(wild_bootstrap_audit, file.path(ROB_OUT, "wild_bootstrap_package_audit.rds"))
saveRDS(wild_bootstrap_results, file.path(ROB_OUT, "wild_bootstrap_results.rds"))
saveRDS(wild_bootstrap_objects, file.path(ROB_OUT, "wild_bootstrap_objects.rds"))
cat("\nWild-cluster-bootstrap audit:\n")
print(wild_bootstrap_audit)
print(wild_bootstrap_results)

# 7. Alternative population denominators (Dominican Republic) ----
pop_main_path <- file.path(DR_CLEAN, "pop_municipio_2016_2025.rds")
pop_b_path <- file.path(DR_CLEAN, "pop_municipio_2016_2025_B.rds")
pop_candidates <- c(series_A = pop_a_path, series_B = pop_b_path, series_C = pop_main_path)

inspect_pop_file <- function(label, path) {
  if (!file.exists(path)) {
    return(data.table(series = label, path = path, exists = FALSE, rows = NA_integer_,
                      localities = NA_integer_, first_year = NA_integer_, last_year = NA_integer_))
  }
  p <- as.data.table(readRDS(path))
  needed <- c("adm3_pcode", "sex", "age_group", "year", "pop")
  if (length(setdiff(needed, names(p))) > 0L) stop("Population file missing required columns: ", path)
  data.table(
    series = label, path = path, exists = TRUE, rows = nrow(p),
    localities = uniqueN(p$adm3_pcode), first_year = min(p$year), last_year = max(p$year)
  )
}

denominator_inventory <- rbindlist(Map(inspect_pop_file, names(pop_candidates), pop_candidates))
denominator_inventory[, interpretation := c(
  "Main: ONE municipal estimates through 2020, extended with Hamilton-Perry and raked to ONE province projections",
  "Robustness: intercensal interpolation between the 2010 and 2022 censuses",
  "Robustness: 2022-census-anchored projection; 158 municipalities before folding to the common 155 geography"
)]
denominator_inventory <- rbind(
  denominator_inventory,
  data.table(
    series = "Ecuador_current", path = ec_path, exists = TRUE, rows = nrow(ec),
    localities = uniqueN(ec$canton), first_year = min(ec$year), last_year = max(ec$year),
    interpretation = "One official-projection denominator set in the final panel; no independent alternative series identified in the repository"
  ),
  fill = TRUE
)
saveRDS(denominator_inventory, file.path(ROB_OUT, "denominator_inventory.rds"))
cat("\nDenominator inventory:\n")
print(denominator_inventory)

build_dr_alt_panel <- function(path, series_label, fold_to_155 = FALSE) {
  p <- as.data.table(readRDS(path))
  if (fold_to_155) {
    p[, adm3_pcode := to_155(adm3_pcode)]
    p <- p[, .(pop = sum(pop)), by = .(adm3_pcode, sex, age_group, year)]
  }
  p <- p[sex == "F" & year %in% dr$year & adm3_pcode %in% dr$id]
  keep_ages <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  p <- p[age_group %in% keep_ages]
  if (p[, anyDuplicated(.SD), .SDcols = c("adm3_pcode", "year", "age_group")] != 0L) {
    stop("Duplicate locality-year-age cells in ", series_label)
  }
  direct <- data.table::dcast(
    p[age_group %in% c("15-19", "20-24", "25-29", "30-34")],
    adm3_pcode + year ~ age_group, value.var = "pop"
  )
  setnames(direct,
           c("adm3_pcode", "15-19", "20-24", "25-29", "30-34"),
           c("id", "pop1519", "pop2024", "pop2529", "pop3034"))
  older <- p[age_group %in% c("35-39", "40-44", "45-49"),
             .(pop3549 = sum(pop)), by = .(id = adm3_pcode, year)]
  den <- merge(direct, older, by = c("id", "year"), all = TRUE)
  base_vars <- c(
    "id", "year", "post", "poorest", XR_DR,
    "n_births_15_19", "n_births_20_24", "n_births_25_29", "n_births_30_34", "nb_3549",
    "womenA_15_19", "womenA_20_24", "womenA_25_29", "womenA_30_34", "women3549"
  )
  d <- merge(dr[, ..base_vars], den, by = c("id", "year"), all.x = TRUE)
  if (nrow(d) != nrow(dr)) stop("Alternative denominator merge changed row count for ", series_label)
  if (d[, anyNA(.SD), .SDcols = paste0("pop", c(AG_TEEN, AG_COMP))]) {
    stop("Alternative denominator merge introduced missing values for ", series_label)
  }
  if (d[, any(unlist(.SD) <= 0), .SDcols = paste0("pop", c(AG_TEEN, AG_COMP))]) {
    stop("Non-positive alternative denominator in ", series_label)
  }

  birth_vars <- c("n_births_15_19", "n_births_20_24", "n_births_25_29", "n_births_30_34", "nb_3549")
  age_codes <- c(AG_TEEN, AG_COMP)
  for (j in seq_along(age_codes)) {
    d[, (paste0("tn", age_codes[j])) := get(birth_vars[j]) / get(paste0("pop", age_codes[j])) * 1000]
    d[, (paste0("log_tn", age_codes[j])) := log(get(paste0("tn", age_codes[j])) + 1)]
  }
  for (k in seq_along(AG_COMP)) {
    d[, (paste0("log_dif", k)) := log_tn1519 - get(paste0("log_tn", AG_COMP[k]))]
  }
  d[, denominator_series := series_label]
  setorder(d, id, year)
  d[]
}

alt_panels <- list(
  series_A = build_dr_alt_panel(pop_a_path, "series_A"),
  series_B = build_dr_alt_panel(pop_b_path, "series_B"),
  series_C = build_dr_alt_panel(pop_main_path, "series_C", fold_to_155 = TRUE)
)
saveRDS(alt_panels, file.path(ROB_OUT, "alternative_denominator_panels.rds"))

alt_denominator_models <- list()
alt_denominator_results <- list()
for (series_name in names(alt_panels)) {
  fitted <- fit_main_models(alt_panels[[series_name]], paste0("DR_", series_name), "id", XR_DR)
  alt_denominator_models <- c(alt_denominator_models, fitted$pooled)
}
alt_denominator_results <- summarize_pooled(alt_denominator_models)
setnames(alt_denominator_results, "country", "series")
alt_denominator_results[, series := sub("^DR_", "", series)]
setorder(alt_denominator_results, series, comparison, specification)
saveRDS(alt_denominator_models, file.path(ROB_OUT, "alternative_denominator_models.rds"))
saveRDS(alt_denominator_results, file.path(ROB_OUT, "alternative_denominator_results.rds"))

baseline_panel <- alt_panels$series_A
den_comparison <- rbindlist(lapply(c("series_B", "series_C"), function(series_name) {
  d <- alt_panels[[series_name]]
  rbindlist(lapply(seq_along(AG_LABS), function(j) {
    alt <- d[[paste0("pop", names(AG_LABS)[j])]]
    baseline <- baseline_panel[[paste0("pop", names(AG_LABS)[j])]]
    pct <- 100 * (alt / baseline - 1)
    data.table(
      series = series_name, age_group = unname(AG_LABS[j]), n = length(pct),
      mean_percent_difference = mean(pct), median_percent_difference = stats::median(pct),
      std_dev_percent_difference = stats::sd(pct), min_percent_difference = min(pct),
      max_percent_difference = max(pct)
    )
  }))
}))
saveRDS(den_comparison, file.path(ROB_OUT, "alternative_denominator_comparison.rds"))
cat("\nAlternative-denominator pooled estimates:\n")
print(alt_denominator_results[, .(series, comparison, specification, estimate, std_error, p_value, n_obs)])
cat("\nAlternative-versus-Series-A denominator differences:\n")
print(den_comparison)

# 8. Results summary ----
pooled_md <- copy(pooled_results)
pooled_md[, `:=`(
  country = country,
  comparison = paste0("15-19 vs ", AG_LABS[comparison]),
  specification = SPEC_LABELS[specification],
  significance = stars(p_value)
)]
setnames(pooled_md,
         c("country", "comparison", "specification", "estimate", "std_error", "p_value", "exact_percent", "n_obs", "significance"),
         c("Country", "Comparison", "Specification", "Estimate", "SE", "p", "Exact percent", "N", "Sig."))

pretrend_md <- copy(pretrend_tests)
pretrend_md[, `:=`(
  comparison = paste0("15-19 vs ", AG_LABS[comparison]),
  specification = SPEC_LABELS[specification]
)]
setnames(pretrend_md,
         c("country", "comparison", "specification", "pre_years", "statistic", "p_value", "df1", "df2", "n_obs"),
         c("Country", "Comparison", "Specification", "Years tested", "F statistic", "p", "df1", "df2", "N"))

placebo_md <- copy(adult_placebo_results)
placebo_md[, `:=`(specification = SPEC_LABELS[specification], significance = stars(p_value))]
setnames(placebo_md,
         c("country", "comparison", "specification", "estimate", "std_error", "p_value", "exact_percent", "n_obs", "significance"),
         c("Country", "Comparison", "Specification", "Estimate", "SE", "p", "Exact percent", "N", "Sig."))

zero_md <- copy(zero_rate_summary)
setnames(zero_md,
         c("country", "age_group", "locality_years", "nonmissing", "missing", "zero_rates", "zero_percent_nonmissing"),
         c("Country", "Age group", "Locality-years", "Nonmissing", "Missing", "Zeros", "Zero percent"))

rate_md <- copy(rate_summary)
setnames(rate_md,
         c("country", "age_group", "n_nonmissing", "n_missing", "mean", "median", "std_dev", "min", "max"),
         c("Country", "Age group", "N", "Missing", "Mean", "Median", "SD", "Min", "Max"))

alt_md <- copy(alt_denominator_results)
alt_md[, `:=`(
  comparison = paste0("15-19 vs ", AG_LABS[comparison]),
  specification = SPEC_LABELS[specification], significance = stars(p_value)
)]
setnames(alt_md,
         c("series", "comparison", "specification", "estimate", "std_error", "p_value", "exact_percent", "n_obs", "significance"),
         c("Series", "Comparison", "Specification", "Estimate", "SE", "p", "Exact percent", "N", "Sig."))

event_md <- copy(event_results)
event_md[, `:=`(
  comparison = paste0("15-19 vs ", AG_LABS[comparison]),
  specification = SPEC_LABELS[specification], significance = stars(p_value)
)]
setnames(event_md,
         c("country", "comparison", "specification", "year", "estimate", "std_error", "p_value", "exact_percent", "n_obs", "significance"),
         c("Country", "Comparison", "Specification", "Year", "Estimate", "SE", "p", "Exact percent", "N", "Sig."))

den_md <- copy(den_comparison)
setnames(den_md,
         c("series", "age_group", "mean_percent_difference", "median_percent_difference", "std_dev_percent_difference", "min_percent_difference", "max_percent_difference", "n"),
         c("Series", "Age group", "Mean pct. diff.", "Median pct. diff.", "SD pct. diff.", "Min pct. diff.", "Max pct. diff.", "N"))

wild_md <- copy(wild_bootstrap_results)
wild_md[, `:=`(
  comparison = paste0("15-19 vs ", AG_LABS[comparison]),
  specification = SPEC_LABELS[specification]
)]
setnames(wild_md,
         c("country", "comparison", "specification", "estimate", "asymptotic_p_value", "bootstrap_p_value", "bootstrap_draws", "significance_changed_5pct"),
         c("Country", "Comparison", "Specification", "Estimate", "Asymptotic p", "WCB p", "Draws", "Changed at 5%"))

summary_lines <- c(
  "# Submission robustness results",
  "",
  paste0("Generated by `analysis/R/17_submission_robustness.R` on ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "."),
  "",
  "## Scope and specification",
  "",
  "The outcome in every main model is the log difference `log(teen rate + 1) - log(comparison rate + 1)`. The pooled model regresses this outcome on Post-by-Poorest, unrestricted year fixed effects, and locality fixed effects. The event study includes the same two fixed effects and year-by-poorest interactions with 2020 omitted. Standard errors are clustered by locality. Models are reported both without covariates and with the registry-derived covariates used in the existing pipeline. Because those covariates may respond to the pandemic, the no-covariate estimates are substantively important rather than merely preliminary.",
  "",
  "`Exact percent` is `100*(exp(beta)-1)`. This is a transformation of the log-difference coefficient; because the underlying rate is transformed as `log(rate+1)`, it is not exactly a conventional percentage change in the untransformed birth rate.",
  "",
  "## Analysis samples",
  "",
  write_md_table(data_inventory, names(data_inventory), digits = 0L),
  "",
  "## Age-specific rate summary statistics",
  "",
  write_md_table(rate_md, c("Country", "Age group", "N", "Missing", "Mean", "Median", "SD", "Min", "Max")),
  "",
  "## Pooled Post x Poorest DDD estimates",
  "",
  write_md_table(pooled_md, c("Country", "Comparison", "Specification", "Estimate", "SE", "p", "Exact percent", "N", "Sig.")),
  "",
  "## Joint Wald tests of all pre-2020 year x Poorest coefficients",
  "",
  "The null is that every listed pre-2020 interaction coefficient is jointly zero. A small p-value is evidence against parallel pretrends for that comparison and specification.",
  "",
  write_md_table(pretrend_md, c("Country", "Comparison", "Specification", "Years tested", "F statistic", "p", "df1", "df2", "N")),
  "",
  "## Adult placebo estimates using the headline log-difference functional form",
  "",
  write_md_table(placebo_md, c("Country", "Comparison", "Specification", "Estimate", "SE", "p", "Exact percent", "N", "Sig.")),
  "",
  "## Zero age-specific locality-year rates",
  "",
  write_md_table(zero_md, c("Country", "Age group", "Locality-years", "Nonmissing", "Missing", "Zeros", "Zero percent")),
  "",
  "`log(rate+1)` keeps zero-rate locality-years in the sample and avoids undefined logarithms. Its cost is that coefficients are not exact proportional changes in the raw rate, especially when rates are close to zero. Sparse zeros are most consequential for comparisons involving age groups with the largest zero share. The RDS file `zero_rate_by_wealth.rds` also reports zeros separately for poorest and other localities.",
  "",
  "## Alternative Dominican population denominators",
  "",
  "The Dominican repository documents three denominator series. Series A is the intended main series: ONE municipal estimates through 2020, extended through 2025 with Hamilton-Perry projections raked to ONE province totals. Series B is an intercensal 2010-2022 series, and Series C is anchored on the 2022 census. Series C is folded from 158 to the common 155-municipality geography here. The saved comparative panel had inadvertently combined Series A for ages 15-34 with an unfolded Series C denominator for ages 35-49; this script reconstructs ages 35-49 on the folded Series A geography before fitting the Dominican main models. No independent Ecuador alternative denominator series was identified.",
  "",
  "### Difference from the Series A main denominators",
  "",
  write_md_table(den_md, c("Series", "Age group", "Mean pct. diff.", "Median pct. diff.", "SD pct. diff.", "Min pct. diff.", "Max pct. diff.", "N")),
  "",
  "### Pooled estimates using each denominator series",
  "",
  write_md_table(alt_md, c("Series", "Comparison", "Specification", "Estimate", "SE", "p", "Exact percent", "N", "Sig.")),
  "",
  "## Wild-cluster-bootstrap inference",
  "",
  paste0("Null-imposed WCR11 tests use Rademacher weights, ", format(WILD_B, big.mark = ","), " bootstrap draws, and the same municipality/canton clustering as the main models. Both base R and `dqrng` seeds are fixed for reproducibility. `fwildclusterboot` version ", utils::packageVersion("fwildclusterboot"), " was used."),
  "",
  write_md_table(wild_md, c("Country", "Comparison", "Specification", "Estimate", "Asymptotic p", "WCB p", "Draws", "Changed at 5%")),
  "",
  "## Full event-study interaction estimates",
  "",
  write_md_table(event_md, c("Country", "Comparison", "Specification", "Year", "Estimate", "SE", "p", "Exact percent", "N", "Sig.")),
  "",
  "## Anomalies and precise limitations",
  "",
  "- The joint Wald tests must be used rather than describing all pre-2020 coefficients as individually small or insignificant. Results that reject the joint null indicate a comparison group with differential pretrends.",
  "- The 20-24 placebo is not null in every Ecuador specification, so ages 20-24 may share part of the post-2020 response and are not a clean unaffected comparison group.",
  "- The oldest comparison group has the most zero locality-year rates. Although `log(rate+1)` retains these cells, it changes the scale and can make estimates sensitive to small counts.",
  "- Registry covariates are constructed from observed births and may be post-treatment variables. Estimates without them avoid that particular conditioning problem but do not solve the design's broader identifying assumptions.",
  "- The pooled model compresses heterogeneous annual effects into one 2021+ average. The event-study estimates should remain the primary description of timing, especially in the Dominican Republic.",
  "- Clustered standard errors use 155 Dominican municipalities and 221 Ecuadorian cantons. The reported WCR11 wild-cluster-bootstrap tests provide a finite-sample inference sensitivity check for the pooled interactions.",
  "- The denominator sensitivity is available only for the Dominican Republic. The Dominican 35-49 outcome in the saved comparative panel mixed denominator series and geographic definitions; this script corrects that outcome to the documented folded Series A baseline before fitting the main Dominican models.",
  "- All models are ecological locality-year analyses. They do not identify individual-level wealth, intended fertility, discrimination, or completed cohort fertility.",
  "",
  "Significance markers: *** p<0.01, ** p<0.05, * p<0.10. Exact estimates, standard errors, p-values, and model objects are preserved as RDS files under `analysis/output/submission_robustness/`."
)
writeLines(summary_lines, SUMMARY_FILE, useBytes = TRUE)

summary_artifacts <- list(
  data_inventory = data_inventory,
  rate_summary = rate_summary,
  pooled_results = pooled_results,
  event_results = event_results,
  pretrend_tests = pretrend_tests,
  adult_placebo_results = adult_placebo_results,
  zero_rate_summary = zero_rate_summary,
  zero_rate_by_wealth = zero_rate_by_wealth,
  wild_bootstrap_audit = wild_bootstrap_audit,
  wild_bootstrap_results = wild_bootstrap_results,
  dr_3549_reconstruction_audit = dr_3549_reconstruction_audit,
  denominator_inventory = denominator_inventory,
  alternative_denominator_results = alt_denominator_results,
  alternative_denominator_comparison = den_comparison
)
saveRDS(summary_artifacts, file.path(ROB_OUT, "results_bundle.rds"))

cat("\nOutput verification:\n")
expected_outputs <- c(
  "data_inventory.rds", "dr_3549_reconstruction_audit.rds", "main_models.rds", "pooled_results.rds", "event_study_results.rds",
  "joint_pretrend_tests.rds", "adult_placebo_models.rds", "adult_placebo_results.rds",
  "rate_summary_statistics.rds", "zero_rate_summary.rds", "zero_rate_by_wealth.rds",
  "wild_bootstrap_package_audit.rds", "wild_bootstrap_results.rds", "wild_bootstrap_objects.rds", "denominator_inventory.rds",
  "alternative_denominator_panels.rds", "alternative_denominator_models.rds",
  "alternative_denominator_results.rds", "alternative_denominator_comparison.rds", "results_bundle.rds"
)
output_paths <- file.path(ROB_OUT, expected_outputs)
print(data.table(file = output_paths, exists = file.exists(output_paths)))
stopifnot(all(file.exists(output_paths)), file.exists(SUMMARY_FILE))

cat("\nCompleted:", format(Sys.time(), tz = "America/Chicago"), "\n")
cat("Results summary:", SUMMARY_FILE, "\n")
cat("RDS directory:", ROB_OUT, "\n")
