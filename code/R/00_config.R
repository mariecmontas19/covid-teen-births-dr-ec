# 00_config.R — shared paths, packages, and plot theme
# Comparative Fertility Paper (DR + Ecuador), MM & OG
# Created 2026-07-27

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(haven)
  library(ggplot2)
})

ROOT    <- Sys.getenv("FERT_ROOT", unset = getwd())  # project folder
DR_CLEAN <- file.path(ROOT, "Dominican Republic", "clean")
DR_CENSUS <- file.path(ROOT, "Dominican Republic", "census")  # 2010 and 2022 census microdata
EC_DATA  <- file.path(ROOT, "Ecuador", "covid_baby_bust_boom", "data")
OUT      <- file.path(ROOT, "analysis", "output")
OUT_TAB  <- file.path(OUT, "tables")
OUT_FIG  <- file.path(OUT, "figures")

# Age-group codes used throughout (mirrors Omar's Ecuador convention)
AG_TEEN  <- "1519"
AG_COMP  <- c("2024", "2529", "3034", "3549")
AG_LABS  <- c("1519" = "15-19", "2024" = "20-24", "2529" = "25-29",
              "3034" = "30-34", "3549" = "35-49")

# Shared figure theme (both countries)
theme_fert <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}
COL_POOR <- "#B2182B"   # poorest quintile
COL_RICH <- "#2166AC"   # rest
