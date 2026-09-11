# MA-ClimateChange-BrainHealth

### Code and data for a meta-analysis of climate change interventions on mental health and well-being outcomes

#### Description

This repository contains the final, harmonized data and analysis code for the Cohen's d meta-analysis of climate-related exposures/interventions and neuropsychological outcomes (overall pool, and separately for mental health and behavioral outcomes).

#### Contents

* **Scripts** (`scripts/`):
  - `domain_utils.R` — shared outcome-domain classification and study-level composite-aggregation helpers, sourced by `run_ma.R`.
  - `run_ma.R` — the meta-analysis engine. Runs the same random-effects Cohen's d analysis once unfiltered ("All") and once per outcome domain (mental_health and behavioral), writing each run's output to its own folder under `results/`.

* **Data** (`data/`):
  - `pool_continuous_d.csv` — the final, harmonized Cohen's d pool actually used by `run_ma.R`. 

* **Software versions used**:
  - Python 3.11.15 — dependencies pinned in `requirements.txt` (pandas, numpy, openpyxl)
  - R 4.3.3 ("Angel Food Cake") — key packages: `meta` 8.2.1, `metafor` 4.8.0, `dmetar` 0.1.0, `tidyverse` 2.0.0, `readxl` 1.4.5, `esc` 0.5.1, `writexl` 1.5.4, `scales` 1.4.0, `patchwork` 1.3.2


#### License

Released under the [MIT License](LICENSE).
