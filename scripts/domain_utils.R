# =============================================================================
# Shared outcome-domain classification & study-level aggregation helpers
# Sourced by Full-MA.R and Mental_Health-MA.R (kept identical between the two
# scripts; previously duplicated verbatim, which let their subgroup/SA
# definitions silently drift apart).
# Requires dplyr (tidyverse) to already be loaded by the calling script.
# =============================================================================

normalize_flag <- function(x) tolower(trimws(as.character(x)))

get_domain <- function(row) {
  d <- c()
  if (isTRUE(row[["cognitive_outcome"]]     == "yes")) d <- c(d, "cognitive")
  if (isTRUE(row[["mental_health_outcome"]] == "yes")) d <- c(d, "mental_health")
  if (isTRUE(row[["brain_outcome"]]         == "yes")) d <- c(d, "brain")
  if (isTRUE(row[["behavioral_outcome"]]    == "yes")) d <- c(d, "behavioral")
  if (length(d) == 0) "unclassified" else paste(d, collapse = "/")
}

# Normalises the four outcome-domain flags, adds outcome_domain (row-level)
# and study_label to pool_d.
prepare_pool_d <- function(pool_d) {
  for (col in c("cognitive_outcome", "mental_health_outcome",
                "brain_outcome", "behavioral_outcome")) {
    pool_d[[col]] <- normalize_flag(pool_d[[col]])
  }
  pool_d$outcome_domain <- apply(pool_d, 1, get_domain)
  pool_d$study_label <- paste0(trimws(pool_d$author), " ", pool_d$year)
  pool_d
}

# =============================================================================
# COMPOSITE EFFECT SIZE (Borenstein et al., 2009)
# =============================================================================
composite_se <- function(se_vec, r = 0.5) {
  se_vec <- se_vec[!is.na(se_vec) & se_vec > 0]
  k <- length(se_vec)
  if (k == 0) return(NA_real_)
  if (k == 1) return(se_vec)
  var_vec   <- se_vec^2
  sum_var   <- sum(var_vec)
  sum_se_sq <- sum(se_vec)^2
  var_c     <- ((1 - r) * sum_var + r * sum_se_sq) / k^2
  sqrt(max(var_c, 0))
}

design_rank <- c("RCT" = 4, "RCT_cluster" = 3,
                 "quasi-experimental" = 2, "observational" = 1)
best_design <- function(designs) {
  ranks <- design_rank[designs]
  ranks[is.na(ranks)] <- 0
  names(which.max(ranks))
}

# Aggregates row-level outcomes (pool_d, already run through prepare_pool_d)
# into one composite effect per study/region/income group. Any group whose
# rows span more than one outcome_domain is labelled "multiple" - its d/se_d
# is a blended composite across domains, not a pure single-domain effect, so
# it must not be treated as belonging to any one domain downstream.
aggregate_d_pool <- function(pool_d) {
  pool_d %>%
    group_by(study_label, region, income_clean) %>%
    summarise(
      k_outcomes       = n(),
      designs_present  = paste(sort(unique(design_family)), collapse = " | "),
      design_family    = best_design(design_family),
      country          = dplyr::first(na.omit(country)),
      outcome_domains  = paste(sort(unique(outcome_domain)), collapse = " | "),
      analysis_fams    = paste(sort(unique(analysis_family)), collapse = " | "),
      # Some studies report a different N per outcome row (different outcomes
      # measured in different sub-samples of the same study). max() is used
      # rather than first() so the reported N reflects the largest measured
      # sub-sample rather than an arbitrary row-order artifact.
      n                = max(N, na.rm = TRUE),
      d                = mean(d_computed, na.rm = TRUE),
      se_d             = composite_se(d_se),
      .groups          = "drop"
    ) %>%
    filter(!is.na(d), !is.na(se_d), se_d > 0) %>%
    mutate(outcome_domain = if_else(grepl("\\|", outcome_domains),
                                    "multiple", outcome_domains))
}
