# COVID-19 and Teen Birth Rates: Comparative Evidence from Ecuador and the Dominican Republic

Replication code for the paper by Marie C. Montás (Harvard T.H. Chan School of Public Health) and Omar Galárraga (Brown University School of Public Health). Contact: mariemontas@fas.harvard.edu.

## Abstract

Although the COVID-19 pandemic produced baby busts across high-income countries, much less is known about how the shock reshaped adolescent fertility in middle-income settings, where teen childbearing is more prevalent and more unequally distributed. Using 3.9 million Ecuadorian and 1.4 million Dominican birth records, we compare adolescent fertility across age groups, time, and localities ranked by household wealth in the 2010 censuses. In Ecuador, conceptions fell 9.4 percent below their projected level between March and December 2020, yet adolescent conceptions in the poorest cantons rose 5.6 percent while those in the remaining 80 percent fell 13.6 percent. In the Dominican Republic, adolescent conceptions initially fell, but by 2023 to 2024 the relative deviations in the poorest municipalities exceeded those elsewhere by 16.5 to 21.1 percentage points. Pooled triple-difference estimates imply that post-2020 adolescent birth rates rose, relative to adult rates and to other localities, by about 13 to 30 percent in the poorest Ecuadorian cantons and 12 to 22 percent in the poorest Dominican municipalities, and Ecuador's adolescent-birth concentration index nearly doubled in magnitude.

## Contents

```
code/R/   analysis pipeline, 00_config.R through 17_submission_robustness.R
```

Data and generated outputs are not included.

## Pipeline

Each script sources `00_config.R`, which takes the project folder from the `FERT_ROOT` environment variable (or the working directory) and defines the data folders relative to it. Run the scripts in numeric order.

| Stage | What it does |
|---|---|
| 00 | Paths, packages, figure theme; 2010-census household asset index for the Dominican Republic (PCA) |
| 01 | Dominican Republic municipality × age × year panel |
| 02-03 | Dominican Republic difference-in-differences, event studies, and preferred specification |
| 04 | Ecuador replication, checked against the coauthor's Stata estimates |
| 05-06 | Figures and tables |
| 07 | HonestDiD sensitivity bounds |
| 08 | Continuous wealth gradient |
| 09-10 | Monthly conceptions and deviations from projected levels, both countries |
| 11 | Paper tables |
| 12 | Sensitivity analyses and placebos |
| 13 | Concentration curves and indices |
| 14 | Mechanisms |
| 15-16 | Household margin and census-based gradient |
| 17 | Robustness checks for submission |

## Software

R 4.3.3 with `data.table` 1.17.8, `fixest` 0.12.1, `haven` 2.5.5, `ggplot2` 3.5.2, `psych` 2.4.3, and `HonestDiD` 0.2.6.

## Data access

Dominican births come from the Ministry of Public Health's live-birth registry (BDNV), provided under a data use agreement, and cannot be redistributed. Ecuadorian births come from INEC's *Estadísticas de Nacidos Vivos*, which are public. Population denominators come from the national statistics offices (ONE and INEC), and the wealth indices are built from the 2010 censuses of both countries.

## Citation

Montás, M. C., and Galárraga, O. (2026). *COVID-19 and Teen Birth Rates: Comparative Evidence from Ecuador and the Dominican Republic.* Working paper.
