# domain_utils.R
# Domain filtering + study-level composite aggregation for the 6 analyses
# (Overall effects + 5 outcome domains). Sourced by run_ma.R, which loads
# dplyr before sourcing this file.
# -------------------------------------------------------------------------------------------

DOMAIN_SPEC <- list(
  overall_effects = list(display = "Overall effects", column = NA_character_),
  brain_health     = list(display = "Brain Health",     column = "Brain Health"),
  mental_health    = list(display = "Mental Health",    column = "Mental Health"),
  general_health   = list(display = "General health",   column = "General health"),
  well_being       = list(display = "Well-being",       column = "Well-being"),
  pro_environmental = list(
    display = "Pro-environmental attitudes, perceptions, and behaviors",
    column  = "Pro-environmental attitudes, perceptions, and behaviors")
)

# Keeps only the rows tagged "yes" for one domain. Overall effects has no
# column to filter on, so it returns the full pool untouched — every
# eligible row once, not a concatenation of the five domains. Domains are
# non-exclusive: a row can be "yes" in more than one column, so the same
# study can show up in more than one domain's filtered pool.
filter_domain <- function(pool, domain_key) {
  spec <- DOMAIN_SPEC[[domain_key]]
  if (is.na(spec$column)) return(pool)
  pool[pool[[spec$column]] == "yes", ]
}

# ── COMPOSITE EFFECT SIZE (Borenstein et al., 2009) ─────────────────────────
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

design_rank <- c("RCT" = 4, "RCT_cluster" = 3, "quasi-experimental" = 2, "observational" = 1)
best_design <- function(designs) {
  ranks <- design_rank[designs]
  ranks[is.na(ranks)] <- 0
  names(which.max(ranks))
}

# Combines every outcome row for the same study/region/income group into one
# composite effect. Composite N is max(N) across that group's rows, not
# sum(N) — a study measuring several outcomes in the same sample shouldn't
# have its population size multiply just because it reported more outcomes.
aggregate_pool <- function(pool) {
  pool %>%
    group_by(study_label, region, income_clean) %>%
    summarise(
      k_outcomes      = n(),
      source_row_ids  = paste(source_row_id, collapse = "|"),
      outcomes_combined = paste(outcome_name, collapse = " // "),
      designs_present = paste(sort(unique(design_family)), collapse = " | "),
      design_family   = best_design(design_family),
      country         = dplyr::first(na.omit(country)),
      n_dz            = sum(is_dz),
      n               = max(N, na.rm = TRUE),
      d               = mean(d_computed, na.rm = TRUE),
      se_d            = composite_se(d_se),
      .groups         = "drop"
    ) %>%
    filter(!is.na(d), !is.na(se_d), se_d > 0)
}

# Same as aggregate_pool(), but excludes within-subject d_z rows before
# building the composites — a sensitivity check for studies whose effect
# size is a pre/post comparison within one group rather than a comparison
# between two groups.
aggregate_pool_no_dz <- function(pool) aggregate_pool(pool[!pool$is_dz, ])

# Same as aggregate_pool(), but keeps only rows where the effect size was
# reported directly in the source article (es_pool == "d_direct"), filtering
# BEFORE building composites rather than after. Filtering first matters: a
# study's composite can blend several outcomes together, so selecting
# "studies with at least one direct row" and then reusing the full composite
# would still let converted effects leak into a supposedly direct-only
# result. Filtering the raw rows first keeps the composite itself pure.
aggregate_pool_direct_only <- function(pool) aggregate_pool(pool[pool$es_pool == "d_direct", ])
