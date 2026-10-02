# =============================================================================
# run_ma.R — standalone meta-analysis pipeline (public repo version)
# =============================================================================
# Loads the harmonized Cohen's d pool (data/pool_continuous_d.csv), then runs
# one parameterized meta-analysis engine once for "Overall effects" and once
# per outcome domain (Brain Health, Mental Health, General health, Well-being,
# Pro-environmental attitudes/perceptions/behaviors — 5 domains, 6 analyses
# total), writing each analysis to its own results/<domain>/ folder, plus a
# joint cross-analysis summary and a full procedure-level run log.
#
# This script is self-contained except for scripts/domain_utils.R (domain
# filtering + study-level composite aggregation), which it sources. Every
# other step — data loading/normalization, the analysis engine, and the
# driver loop — lives in this one file, so the public repo doesn't need the
# multi-file code/ layout the working pipeline uses internally.
#
# Every procedure (forest, outliers, funnel/egger, GOSH, subgroups, meta-
# regression, sensitivity analyses) is isolated via safe_run(): an error or
# warning in one procedure is logged and does NOT stop the remaining
# procedures, the remaining moderators, or the remaining domains.
#
# Reproducible entry point:
#   cd MA-ClimateChange-BrainHealth && Rscript scripts/run_ma.R
# (also runs correctly via `Rscript scripts/run_ma.R` from any directory —
# the project root is resolved from this script's own path.)
# =============================================================================
rm(list = ls()); gc()

suppressMessages({
  library(dplyr)
  library(meta)
  library(metafor)
  library(dmetar)
  library(writexl)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
PROJECT_ROOT <- if (length(script_arg) > 0) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
} else {
  normalizePath(getwd())
}

source(file.path(PROJECT_ROOT, "scripts", "domain_utils.R"))

# =============================================================================
# STAGE 1 — load + validate the pool (data/pool_continuous_d.csv)
# =============================================================================
SOURCE_CSV <- file.path(PROJECT_ROOT, "data", "pool_continuous_d.csv")

if (!file.exists(SOURCE_CSV)) {
  stop("Can't find ", SOURCE_CSV, ". Run this from the project root, e.g.:\n",
       "  cd MA-ClimateChange-BrainHealth && Rscript scripts/run_ma.R")
}

DOMAIN_COLS <- c(
  "Brain Health", "Mental Health", "General health", "Well-being",
  "Pro-environmental attitudes, perceptions, and behaviors"
)

normalize_region <- function(x) {
  x <- trimws(as.character(x))
  x[x == "north America"] <- "North America"
  x
}

normalize_income <- function(x) {
  x <- tolower(trimws(x))
  case_when(
    x == "high income"         ~ "High income",
    x == "upper-middle income" ~ "Upper-middle income",
    x == "lower-middle income" ~ "Lower-middle income",
    TRUE                       ~ x  # unexpected value: kept as-is, not merged into a category
  )
}

clean_year <- function(y) {
  y <- trimws(as.character(y))
  ifelse(grepl("^[0-9]+\\.0$", y), sub("\\.0$", "", y), y)  # "2025.0"->"2025"; "2013b" untouched
}

normalize_domain_flag <- function(x) tolower(trimws(as.character(x)))

is_within_subject_dz <- function(effect_size_type_orig) {
  grepl("d_z|within", effect_size_type_orig, ignore.case = TRUE)
}

# Reads the harmonized pool and normalizes it in memory (region/income
# casing, year suffix, domain-flag text, numeric types). Does not recompute
# effect sizes or standard errors; the pool already arrives with those
# computed.
load_and_validate_pool <- function(csv_path = SOURCE_CSV) {
  stopifnot(file.exists(csv_path))

  # colClasses="character" on first pass so nothing silently coerces (e.g.
  # "2013b" would break a numeric column read); numeric columns are cast
  # explicitly below. check.names=FALSE preserves domain column names that
  # contain spaces/commas (e.g. "Brain Health") as-is.
  df <- read.csv(csv_path, colClasses = "character", check.names = FALSE,
                  stringsAsFactors = FALSE, na.strings = "", encoding = "UTF-8")

  required_cols <- c("author", "year", "N", "n_exp", "n_ctrl", "analysis_n",
                      "design_family", "region", "country", "income_status",
                      "outcome_name", "outcome_type", DOMAIN_COLS,
                      "predictor_category", "d_computed", "d_se",
                      "transformed_effect_type", "original_effect_size",
                      "effect_size_type_orig", "es_pool", "conversion_method",
                      "effect_size_SE", "QC_flag", "analysis_family")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df$source_row_id <- sprintf("row_%04d", seq_len(nrow(df)))

  # ── Numeric casts ────────────────────────────────────────────────────────
  num_cols <- c("N", "n_exp", "n_ctrl", "analysis_n", "d_computed", "d_se",
                "original_effect_size", "effect_size_SE")
  for (col in num_cols) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))

  # ── Text/locale normalization (in this derived copy only) ───────────────
  df$region <- normalize_region(df$region)
  df$income_clean <- normalize_income(df$income_status)

  df$year_clean <- clean_year(df$year)
  df$study_label <- paste(trimws(df$author), df$year_clean)

  for (col in DOMAIN_COLS) df[[col]] <- normalize_domain_flag(df[[col]])

  df$design_family <- trimws(df$design_family)

  # ── Derived fields needed by downstream aggregation/sensitivity analyses ──
  df$is_dz <- is_within_subject_dz(df$effect_size_type_orig)

  df
}

save_derived_pool <- function(df, out_dir = file.path(PROJECT_ROOT, "data_derived")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(df, file.path(out_dir, "pool_prepared.csv"), row.names = FALSE, na = "")
}

# =============================================================================
# STAGE 2 — the analysis engine (run once per analysis: Overall effects + 5
# domains)
# =============================================================================

# ── Run log: one row per attempted procedure, across all domains ───────────
RUN_LOG <- new.env()
RUN_LOG$rows <- list()

log_row <- function(domain, analysis, moderator = NA_character_, filter_used = NA_character_,
                     k_used = NA_integer_, n_used = NA_integer_, status, message = NA_character_,
                     output_path = NA_character_) {
  RUN_LOG$rows[[length(RUN_LOG$rows) + 1]] <- data.frame(
    domain = domain, analysis = analysis, moderator = moderator, filter_used = filter_used,
    k_used = as.integer(k_used), n_used = as.integer(n_used), status = status,
    message = ifelse(is.na(message), "", message),
    output_path = ifelse(is.na(output_path), "", output_path),
    timestamp = format(Sys.time()),
    stringsAsFactors = FALSE
  )
}

# Executes expr_fun() with error+warning isolation. Guarantees any open sink()
# or graphics device gets closed even if expr_fun errors partway through a
# plot, so one failed figure can't corrupt subsequent sink()/plot calls.
safe_run <- function(domain, analysis, expr_fun, moderator = NA_character_, filter_used = NA_character_,
                      k_used = NA_integer_, n_used = NA_integer_, output_path = NA_character_) {
  status <- "completed"; msg <- NA_character_; warn_msgs <- character(0); result <- NULL
  result <- tryCatch({
    withCallingHandlers(
      expr_fun(),
      warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
  }, error = function(e) {
    status <<- "error"; msg <<- conditionMessage(e); NULL
  })
  while (sink.number() > 0) sink()
  while (dev.cur() != 1) dev.off()
  if (identical(status, "completed") && length(warn_msgs) > 0) {
    status <- "completed_with_warnings"
    msg <- paste(unique(warn_msgs), collapse = " | ")
  }
  log_row(domain, analysis, moderator, filter_used, k_used, n_used, status, msg, output_path)
  result
}

log_not_estimable <- function(domain, analysis, moderator = NA_character_, filter_used = NA_character_,
                               k_used = NA_integer_, message) {
  log_row(domain, analysis, moderator, filter_used, k_used, NA_integer_, "not_estimable", message, NA_character_)
}
log_blocked <- function(domain, analysis, moderator = NA_character_, message) {
  log_row(domain, analysis, moderator, status = "blocked", message = message)
}

sa_val <- function(obj, field) { if (is.null(obj)) return(NA_real_); v <- obj[[field]]; if (is.null(v) || length(v) == 0) NA_real_ else v }
# meta::metagen stores I2 as a fraction (0-1), not a percentage — the printed
# summary multiplies by 100 for display ("I^2 = 94.8%"). Any output table
# that reports I2 alongside the printed .txt files must do the same
# conversion, or it silently shows e.g. 0.9 instead of 94.8.
sa_val_i2_pct <- function(obj) { v <- sa_val(obj, "I2"); if (is.na(v)) NA_real_ else v * 100 }
sa_k   <- function(obj) { if (is.null(obj)) return(NA_integer_); k <- obj[["k"]]; if (is.null(k) || length(k) == 0) NA_integer_ else k }

# ── Output-writing helpers (no internal tryCatch — safe_run() at the call
# site provides isolation + device/sink cleanup) ─────────────────────────────
save_txt <- function(save_path, obj, filename) {
  sink(file.path(save_path, filename)); print(obj); sink()
}

make_forest <- function(save_path, model, filename_pdf, xlim, xlab) {
  # cairo_pdf (not base pdf()) — base pdf()'s default font encoding can't
  # render some study labels' non-ASCII characters and silently substitutes
  # a dot for each one instead. Height scales with k (same heuristic as
  # run_subgroup()'s plot_height) — a fixed page height clips the header row
  # and the bottom statistics block once a domain has enough studies.
  plot_height <- max(9, 0.3 * model$k + 4)
  cairo_pdf(file.path(save_path, filename_pdf), width = 11, height = plot_height, bg = "white")
  forest(model, sortvar = model$TE, smlab = "", leftcols = c("studlab"),
         prediction = TRUE, zero.pval = TRUE, xlim = xlim, xlab = xlab,
         text.random = "Random effects model", text.predict = "Prediction interval",
         addrows.below.overall = 2, print.tau2 = TRUE, print.I2 = TRUE,
         test.overall = TRUE, rightcols = c("ci", "w.random"),
         rightlabs = c("95% CI", "Weight"), col.square = "#2B7CB5",
         col.diamond = "#D94F3D", fontsize = 10)
  dev.off()
}

make_funnel <- function(save_path, model, filename_pdf) {
  col.contour <- c("gray75", "gray85", "gray95")
  pdf(file.path(save_path, filename_pdf), width = 8, height = 7, bg = "white")
  funnel(model, yaxs = "i", lwd = 1, cex = 1.2, contour = c(0.9, 0.95, 0.99), col.contour = col.contour)
  legend("topright", legend = c("p < 0.1", "p < 0.05", "p < 0.01"), fill = col.contour, bty = "n")
  dev.off()
}

# GOSH plots visualize heterogeneity across every possible study subset.
# Saved as PNG rather than a vector PDF — with a few dozen composites the
# plot has thousands of points, which renders as tens of MB in PDF but stays
# small as PNG. Set SKIP_GOSH=1 to skip this step (it's the slowest one).
run_gosh <- function(save_path, model, filename_png) {
  rma_obj <- rma(yi = model$TE, sei = model$seTE, method = model$method.tau, test = "knha")
  res_gosh <- gosh(rma_obj, parallel = "multicore", ncpus = 2, progbar = FALSE)
  xr <- range(model$TE, na.rm = TRUE)
  xlim <- c(min(-2, xr[1] - 0.5), max(2, xr[2] + 0.5))
  png(file.path(save_path, filename_png), width = 8, height = 7, units = "in", res = 300, bg = "white")
  plot(res_gosh, alpha = 0.4, het = "I2", col = c("gray70", "#2B7CB5"), xlim = xlim)
  dev.off()
  res_gosh
}

run_subgroup <- function(save_path, base_model, subgroup_var, label, prefix, format = "pdf") {
  model_sg <- update(base_model, subgroup = subgroup_var)
  save_txt(save_path, model_sg, paste0(prefix, "_subgroup_", label, ".txt"))
  plot_file <- file.path(save_path, paste0(prefix, "_subgroup_", label, ".", format))
  plot_height <- max(8, 0.3 * nrow(base_model$data) + 4)
  if (format == "png") {
    png(plot_file, width = 13, height = plot_height, units = "in", res = 300, bg = "white")
  } else {
    pdf(plot_file, width = 13, height = plot_height, bg = "white")
  }
  forest(model_sg, sortvar = model_sg$TE, prediction = TRUE,
         print.tau2 = TRUE, print.I2 = TRUE, test.subgroup = TRUE,
         leftcols = c("studlab"), rightcols = c("ci", "w.random"), rightlabs = c("95% CI", "Weight"),
         col.square = "#2B7CB5", col.diamond = "#D94F3D", fontsize = 9)
  dev.off()
  model_sg
}

run_metareg <- function(save_path, data, formula_str, label, prefix, bubble_mod = NULL) {
  rma_obj <- rma(yi = data$d, sei = data$se_d, mods = as.formula(formula_str),
                  method = "REML", test = "knha", data = data)
  sink(file.path(save_path, paste0(prefix, "_metareg_", label, ".txt")))
  cat("Meta-regression:", label, "\nFormula:", formula_str, "\n\n")
  print(summary(rma_obj))
  cat("\nR^2:", round(rma_obj$R2, 2), "%\n")
  sink()
  if (!is.null(bubble_mod)) {
    pdf(file.path(save_path, paste0(prefix, "_bubble_", label, ".pdf")), width = 7, height = 6, bg = "white")
    regplot(rma_obj, mod = bubble_mod, xlab = label, ylab = "Effect size", col = "#2B7CB5", bg = "#A8C8E8", pch = 21)
    dev.off()
  }
  rma_obj
}

# Releveling to a substantively meaningful reference category (e.g.
# "observational" for design, "Europe" for region) rather than R's default
# alphabetical first level, so the meta-regression coefficients read as a
# comparison against the category readers actually expect.
run_moderator_metareg <- function(save_path, pool_agg, moderator_col, label, prefix, preferred_reference = NULL) {
  levels_present <- levels(factor(pool_agg[[moderator_col]]))
  ref <- if (!is.null(preferred_reference) && preferred_reference %in% levels_present) preferred_reference else levels_present[1]
  fct_col <- paste0(moderator_col, "_fct")
  pool_agg[[fct_col]] <- relevel(factor(pool_agg[[moderator_col]]), ref = ref)
  run_metareg(save_path, pool_agg, paste0("~ ", fct_col), label, prefix)
}

# ── Full pipeline for ONE analysis (Overall effects, or one of the 5 domains) ──
run_domain_ma <- function(pool_full, domain_key, save_path, region_dir, country_dir) {
  label <- DOMAIN_SPEC[[domain_key]]$display
  # Wipe any previous run's output for this domain first — otherwise a file
  # written by an earlier run can survive untouched if the same analysis
  # becomes not_estimable on a later run (a moderator category disappearing,
  # etc.), and get silently picked up as stale data downstream.
  unlink(save_path, recursive = TRUE)
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

  pool <- filter_domain(pool_full, domain_key)
  log_row(label, "row_filter", k_used = nrow(pool), status = "completed",
          message = sprintf("%d source rows eligible for this domain.", nrow(pool)))

  pool_agg <- aggregate_pool(pool)
  k <- nrow(pool_agg)
  write.csv(pool_agg, file.path(save_path, "00_study_summary.csv"), row.names = FALSE)

  if (k < 2) {
    log_not_estimable(label, "base_model", k_used = k,
                       message = sprintf("Only %d independent study composite(s) after aggregation — random-effects model needs >=2. See 00_study_summary.csv.", k))
    for (a in c("forest", "outliers", "funnel_egger", "pcurve", "gosh",
                "subgroup_design", "subgroup_region", "subgroup_income",
                "metareg_design", "metareg_region", "metareg_income", "metareg_log_n",
                "leave_one_out", "sensitivity_rct", "sensitivity_observational",
                "sensitivity_direct_only", "sensitivity_no_dz", "tau_estimator_comparison")) {
      log_blocked(label, a, message = "Blocked: base model not estimable (k < 2).")
    }
    return(invisible(NULL))
  }
  if (k < 10) {
    log_row(label, "small_k_warning", k_used = k, status = "completed_with_warnings",
            message = sprintf("k = %d < 10: interpret heterogeneity/subgroup/moderator results with caution. Every estimable procedure still runs — small k is a warning, not a block.", k))
  }

  n_by_region <- pool_agg %>% group_by(region) %>% summarise(k = n(), N = sum(n, na.rm = TRUE), .groups = "drop")
  write.csv(n_by_region, file.path(region_dir, paste0("N_by_region_", domain_key, ".csv")), row.names = FALSE)
  n_by_country <- pool_agg %>% group_by(country) %>% summarise(k = n(), N = sum(n, na.rm = TRUE), .groups = "drop")
  write.csv(n_by_country, file.path(country_dir, paste0("N_by_country_", domain_key, ".csv")), row.names = FALSE)

  # ── 1. Base model ──────────────────────────────────────────────────────────
  meta_model <- safe_run(label, "base_model", k_used = k, n_used = sum(pool_agg$n, na.rm = TRUE),
    output_path = file.path(save_path, "01_meta_base.txt"),
    expr_fun = function() {
      m <- metagen(TE = d, seTE = se_d, studlab = study_label, data = pool_agg,
                    sm = "SMD", method.tau = "REML", prediction = TRUE,
                    common = FALSE, random = TRUE, method.random.ci = "HK",
                    title = paste("Meta-analysis (SMD) -", label))
      save_txt(save_path, m, "01_meta_base.txt")
      m
    })
  if (is.null(meta_model)) {
    log_blocked(label, "downstream_all", message = "Base model failed to fit (see base_model row) — all downstream procedures blocked for this domain.")
    return(invisible(NULL))
  }

  safe_run(label, "forest", k_used = k, output_path = file.path(save_path, "02_forest.pdf"),
           expr_fun = function() make_forest(save_path, meta_model, "02_forest.pdf", c(-3, 3), "Cohen's d"))

  # ── 2. Outliers ─────────────────────────────────────────────────────────────
  outliers <- safe_run(label, "outliers", k_used = k, output_path = file.path(save_path, "03_outliers.txt"),
    expr_fun = function() { o <- find.outliers(meta_model); save_txt(save_path, o, "03_outliers.txt"); o })
  meta_no_out <- meta_model
  if (!is.null(outliers)) {
    meta_no_out <- safe_run(label, "no_outliers_model", k_used = k, output_path = file.path(save_path, "03_meta_no_outliers.txt"),
      expr_fun = function() {
        m <- update(meta_model, exclude = meta_model$studlab %in% outliers$out.study.random)
        save_txt(save_path, m, "03_meta_no_outliers.txt")
        make_forest(save_path, m, "03_forest_no_outliers.pdf", c(-3, 3), "Cohen's d")
        m
      })
    if (is.null(meta_no_out)) meta_no_out <- meta_model
  } else {
    log_not_estimable(label, "no_outliers_model", k_used = k, message = "Outlier detection did not return a usable result; no-outliers model skipped.")
  }

  # ── 3. Publication bias ──────────────────────────────────────────────────────
  safe_run(label, "funnel", k_used = k, output_path = file.path(save_path, "04_funnel.pdf"),
           expr_fun = function() make_funnel(save_path, meta_model, "04_funnel.pdf"))
  safe_run(label, "egger", k_used = k, output_path = file.path(save_path, "04_egger.txt"),
           expr_fun = function() save_txt(save_path, eggers.test(meta_model), "04_egger.txt"))
  if (!identical(meta_no_out, meta_model)) {
    safe_run(label, "egger_no_outliers", k_used = k, output_path = file.path(save_path, "04_egger_no_out.txt"),
             expr_fun = function() save_txt(save_path, eggers.test(meta_no_out), "04_egger_no_out.txt"))
  }
  safe_run(label, "pcurve", k_used = k, output_path = file.path(save_path, "04_pcurve.pdf"),
    expr_fun = function() {
      pdf(file.path(save_path, "04_pcurve.pdf"), width = 7, height = 7, bg = "white")
      sink(file.path(save_path, "04_pcurve.txt")); print(pcurve(meta_model)); sink()
      dev.off()
    })

  # ── 4. GOSH — always attempted, PNG only (skip if SKIP_GOSH=1 env var set) ───
  if (Sys.getenv("SKIP_GOSH", "0") == "1") {
    log_row(label, "gosh", k_used = k, status = "completed_with_warnings",
            message = "Skipped: SKIP_GOSH=1 (compute-expensive, disabled for this run).")
  } else {
    safe_run(label, "gosh", k_used = k, output_path = file.path(save_path, "05_gosh.png"),
             expr_fun = function() run_gosh(save_path, meta_model, "05_gosh.png"))
  }

  # ── 5. Subgroups (design/region/income) ──────────────────────────────────────
  run_subgroup_checked <- function(colname, sublabel, format = "pdf") {
    vals <- pool_agg[[colname]]
    present <- unique(na.omit(vals))
    if (length(present) < 2) {
      log_not_estimable(label, paste0("subgroup_", sublabel), moderator = colname, k_used = k,
                         message = sprintf("Only %d level(s) present (%s) — subgroup comparison needs >=2.", length(present), paste(present, collapse = ", ")))
      return(invisible(NULL))
    }
    safe_run(label, paste0("subgroup_", sublabel), moderator = colname, k_used = k,
             output_path = file.path(save_path, paste0("06_subgroup_", sublabel, ".txt")),
             expr_fun = function() run_subgroup(save_path, meta_model, vals, sublabel, "06", format = format))
  }
  run_subgroup_checked("design_family", "design", format = "png")  # design subgroup stays PNG; region/income are PDF
  run_subgroup_checked("region", "region")
  run_subgroup_checked("income_clean", "income")

  # ── 6. Meta-regression (design/region/income/log(N)) ─────────────────────────
  run_metareg_checked <- function(colname, sublabel, ref = NULL) {
    levels_present <- unique(na.omit(pool_agg[[colname]]))
    if (length(levels_present) < 2) {
      log_not_estimable(label, paste0("metareg_", sublabel), moderator = colname, k_used = k,
                         message = sprintf("Only %d level(s) present — meta-regression needs >=2.", length(levels_present)))
      return(invisible(NULL))
    }
    safe_run(label, paste0("metareg_", sublabel), moderator = colname, k_used = k,
             output_path = file.path(save_path, paste0("07_metareg_", sublabel, ".txt")),
             expr_fun = function() run_moderator_metareg(save_path, pool_agg, colname, sublabel, "07", preferred_reference = ref))
  }
  run_metareg_checked("design_family", "design", ref = if ("observational" %in% pool_agg$design_family) "observational" else NULL)
  run_metareg_checked("region", "region", ref = if ("Europe" %in% pool_agg$region) "Europe" else NULL)
  run_metareg_checked("income_clean", "income", ref = if ("High income" %in% pool_agg$income_clean) "High income" else NULL)

  if (n_distinct(pool_agg$n) < 2) {
    log_not_estimable(label, "metareg_log_n", moderator = "log(n)", k_used = k,
                       message = "n is constant across composites — log(N) meta-regression needs variation.")
  } else {
    safe_run(label, "metareg_log_n", moderator = "log(n)", k_used = k,
             output_path = file.path(save_path, "07_metareg_log_n.txt"),
             expr_fun = function() run_metareg(save_path, pool_agg, "~ log(n)", "log_n", "07", bubble_mod = 2))
  }

  # Comparing the six analyses against each other as if they were subgroups
  # of Overall effects isn't reproduced: the five domain columns are
  # non-exclusive, so the same composite can carry "yes" in more than one
  # domain at once. A subgroup/omnibus test assumes each composite sits in
  # exactly one bucket — forcing that here would mean either dropping
  # multi-domain composites or silently mislabeling them as single-domain,
  # both of which would misrepresent the comparison. Logged as not_estimable
  # with the reason recorded, rather than forced.
  if (identical(domain_key, "overall_effects")) {
    log_not_estimable(label, "cross_domain_comparison", k_used = k,
      message = "Not reproduced: the 5 domain columns are non-exclusive (a composite can be 'yes' in >1 domain simultaneously), so a single-bucket subgroup/omnibus comparison across domains would violate the independent-groups assumption of meta::subgroup().")
  }

  # ── 7. Leave-one-out ──────────────────────────────────────────────────────────
  safe_run(label, "leave_one_out", k_used = k, output_path = file.path(save_path, "08_loo.pdf"),
    expr_fun = function() {
      loo <- metainf(meta_model)
      loo_height <- max(8, 0.3 * k + 4)
      cairo_pdf(file.path(save_path, "08_loo.pdf"), width = 10, height = loo_height, bg = "white")
      forest(loo, xlab = "Cohen's d", leftcols = c("studlab"), col.square = "#2B7CB5", fontsize = 9)
      dev.off()
    })

  # ── 8. Sensitivity analyses ────────────────────────────────────────────────
  design_present <- unique(na.omit(pool_agg$design_family))
  meta_rct <- NULL
  if (any(c("RCT", "RCT_cluster") %in% design_present)) {
    meta_rct <- safe_run(label, "sensitivity_rct", k_used = sum(pool_agg$design_family %in% c("RCT", "RCT_cluster")),
      output_path = file.path(save_path, "09_sa_rct_only.txt"),
      expr_fun = function() {
        m <- update(meta_model, subset = pool_agg$design_family %in% c("RCT", "RCT_cluster"))
        save_txt(save_path, m, "09_sa_rct_only.txt"); m
      })
  } else {
    log_not_estimable(label, "sensitivity_rct", k_used = 0, message = "No RCT/RCT_cluster composites present in this domain.")
  }

  meta_obs <- NULL
  if ("observational" %in% design_present) {
    meta_obs <- safe_run(label, "sensitivity_observational", k_used = sum(pool_agg$design_family == "observational"),
      output_path = file.path(save_path, "09_sa_observational_only.txt"),
      expr_fun = function() {
        m <- update(meta_model, subset = pool_agg$design_family == "observational")
        save_txt(save_path, m, "09_sa_observational_only.txt"); m
      })
  } else {
    log_not_estimable(label, "sensitivity_observational", k_used = 0, message = "No observational composites present in this domain.")
  }

  # Direct-effects-only: filters raw rows to es_pool=="d_direct" before
  # aggregating, then builds a fresh model, rather than subsetting the
  # already-built composites — a composite can blend several outcomes
  # together, so subsetting after the fact could still leave converted
  # effects mixed into a supposedly direct-only result.
  n_direct_rows <- sum(pool$es_pool == "d_direct")
  meta_direct <- NULL
  if (n_direct_rows > 0) {
    meta_direct <- safe_run(label, "sensitivity_direct_only", filter_used = "es_pool == d_direct, filtered before aggregation",
      k_used = n_direct_rows, output_path = file.path(save_path, "09_sa_direct_only.txt"),
      expr_fun = function() {
        agg_direct <- aggregate_pool_direct_only(pool)
        if (nrow(agg_direct) < 2) stop(sprintf("Only %d composite(s) built from d_direct-only rows — needs >=2.", nrow(agg_direct)))
        m <- metagen(TE = d, seTE = se_d, studlab = study_label, data = agg_direct,
                      sm = "SMD", method.tau = "REML", prediction = TRUE,
                      common = FALSE, random = TRUE, method.random.ci = "HK",
                      title = paste("Direct effects only -", label))
        save_txt(save_path, m, "09_sa_direct_only.txt"); m
      })
  } else {
    log_not_estimable(label, "sensitivity_direct_only", k_used = 0, message = "No rows with es_pool == 'd_direct' in this domain.")
  }

  # Extra sensitivity: excludes within-subject d_z-type rows before
  # aggregating, to check whether results hold when only between-group
  # comparisons are used.
  n_nondz_rows <- sum(!pool$is_dz)
  meta_no_dz <- NULL
  if (sum(pool$is_dz) > 0 && n_nondz_rows > 0) {
    meta_no_dz <- safe_run(label, "sensitivity_no_dz", filter_used = "is_dz == FALSE, filtered before aggregation",
      k_used = n_nondz_rows, output_path = file.path(save_path, "09_sa_no_dz.txt"),
      expr_fun = function() {
        agg_nodz <- aggregate_pool_no_dz(pool)
        if (nrow(agg_nodz) < 2) stop(sprintf("Only %d composite(s) after excluding d_z rows — needs >=2.", nrow(agg_nodz)))
        m <- metagen(TE = d, seTE = se_d, studlab = study_label, data = agg_nodz,
                      sm = "SMD", method.tau = "REML", prediction = TRUE,
                      common = FALSE, random = TRUE, method.random.ci = "HK",
                      title = paste("Excl. within-subject d_z -", label))
        save_txt(save_path, m, "09_sa_no_dz.txt"); m
      })
  } else {
    log_not_estimable(label, "sensitivity_no_dz", k_used = 0, message = "No within-subject d_z rows in this domain — sensitivity is identical to the primary model, not run separately.")
  }

  # Tau-estimator comparison — always attempted, regardless of which
  # estimator is primary. If estimators diverge meaningfully, that's a
  # finding worth reporting, not something to pick around.
  tau_comparison <- safe_run(label, "tau_estimator_comparison", k_used = k,
    output_path = file.path(save_path, "09_sa_tau_estimators.csv"),
    expr_fun = function() {
      comp <- data.frame(estimator = c("PM", "REML", "DL", "HS"), tau2 = NA_real_, I2 = NA_real_, TE = NA_real_)
      for (i in seq_along(comp$estimator)) {
        m <- tryCatch(update(meta_model, method.tau = comp$estimator[i], method.random.ci = "classic"), error = function(e) NULL)
        if (!is.null(m)) { comp$tau2[i] <- m$tau2; comp$I2[i] <- m$I2 * 100; comp$TE[i] <- m$TE.random }
      }
      write.csv(comp, file.path(save_path, "09_sa_tau_estimators.csv"), row.names = FALSE)
      comp
    })

  # ── 9. Sensitivity summary table + comparison forest (PNG) ──────────────────
  safe_run(label, "sensitivity_summary", k_used = k, output_path = file.path(save_path, "09_sensitivity_summary.csv"),
    expr_fun = function() {
      sa_summary <- data.frame(
        analysis = c("Base (all composites)", "Without outliers", "RCT/RCT_cluster only", "Observational only", "Direct d only", "Excl. within-subject d_z"),
        k        = c(sa_k(meta_model), sa_k(meta_no_out), sa_k(meta_rct), sa_k(meta_obs), sa_k(meta_direct), sa_k(meta_no_dz)),
        d        = round(c(sa_val(meta_model,"TE.random"), sa_val(meta_no_out,"TE.random"), sa_val(meta_rct,"TE.random"), sa_val(meta_obs,"TE.random"), sa_val(meta_direct,"TE.random"), sa_val(meta_no_dz,"TE.random")), 3),
        ci_lower = round(c(sa_val(meta_model,"lower.random"), sa_val(meta_no_out,"lower.random"), sa_val(meta_rct,"lower.random"), sa_val(meta_obs,"lower.random"), sa_val(meta_direct,"lower.random"), sa_val(meta_no_dz,"lower.random")), 3),
        ci_upper = round(c(sa_val(meta_model,"upper.random"), sa_val(meta_no_out,"upper.random"), sa_val(meta_rct,"upper.random"), sa_val(meta_obs,"upper.random"), sa_val(meta_direct,"upper.random"), sa_val(meta_no_dz,"upper.random")), 3),
        I2       = round(c(sa_val_i2_pct(meta_model), sa_val_i2_pct(meta_no_out), sa_val_i2_pct(meta_rct), sa_val_i2_pct(meta_obs), sa_val_i2_pct(meta_direct), sa_val_i2_pct(meta_no_dz)), 1),
        p_value  = round(c(sa_val(meta_model,"pval.random"), sa_val(meta_no_out,"pval.random"), sa_val(meta_rct,"pval.random"), sa_val(meta_obs,"pval.random"), sa_val(meta_direct,"pval.random"), sa_val(meta_no_dz,"pval.random")), 4)
      )
      write.csv(sa_summary, file.path(save_path, "09_sensitivity_summary.csv"), row.names = FALSE)

      png(file.path(save_path, "09_sa_comparison_forest.png"), width = 9, height = 6, units = "in", res = 300, bg = "white")
      par(mar = c(5, max(nchar(sa_summary$analysis)) * 0.6 + 1, 4, 1))
      valid <- !is.na(sa_summary$d)
      plot(x = sa_summary$d[valid], y = (nrow(sa_summary):1)[valid],
           xlim = c(min(sa_summary$ci_lower, na.rm = TRUE) - 0.1, max(sa_summary$ci_upper, na.rm = TRUE) + 0.1),
           yaxt = "n", xlab = "Cohen's d (RE)", ylab = "", pch = 18, cex = 1.5, col = "#2B7CB5", main = "")
      segments(x0 = sa_summary$ci_lower[valid], x1 = sa_summary$ci_upper[valid],
                y0 = (nrow(sa_summary):1)[valid], y1 = (nrow(sa_summary):1)[valid], col = "#2B7CB5", lwd = 2)
      axis(2, at = nrow(sa_summary):1, labels = sa_summary$analysis, las = 1, cex.axis = 0.8)
      abline(v = 0, lty = 2, col = "gray50")
      dev.off()
      sa_summary
    })

  invisible(meta_model)
}

# ── Joint no-outliers summary — one row per analysis, built from each
# domain's outliers-excluded model (dmetar::find.outliers()). Refits the
# base + outlier-exclusion models directly (fast — no GOSH, no plots)
# instead of threading meta_no_out out of run_domain_ma()'s return value, so
# it can be called on its own after the main loop.
build_no_outliers_summary <- function(pool_full) {
  rows <- list()
  for (dk in names(DOMAIN_SPEC)) {
    label <- DOMAIN_SPEC[[dk]]$display
    pool <- filter_domain(pool_full, dk)
    agg <- aggregate_pool(pool)
    k <- nrow(agg)
    if (k < 2) {
      rows[[dk]] <- data.frame(analysis = label, k = NA_integer_, N = NA_integer_, d = NA_real_,
                 ci_lower = NA_real_, ci_upper = NA_real_, I2 = NA_real_, tau2 = NA_real_,
                 p_value = NA_real_, n_outliers_excluded = NA_integer_, status = "not_estimable")
      next
    }
    meta_model <- metagen(TE = d, seTE = se_d, studlab = study_label, data = agg, sm = "SMD",
                           method.tau = "REML", prediction = TRUE, common = FALSE, random = TRUE,
                           method.random.ci = "HK", title = label)
    outliers <- tryCatch(find.outliers(meta_model), error = function(e) NULL)
    meta_no_out <- meta_model
    n_excl <- 0L
    if (!is.null(outliers) && length(outliers$out.study.random) > 0) {
      meta_no_out <- tryCatch(update(meta_model, exclude = meta_model$studlab %in% outliers$out.study.random),
                               error = function(e) meta_model)
      n_excl <- length(outliers$out.study.random)
    }
    n_used <- if (n_excl > 0) sum(agg$n[!(agg$study_label %in% outliers$out.study.random)], na.rm = TRUE) else sum(agg$n, na.rm = TRUE)
    rows[[dk]] <- data.frame(
      analysis = label, k = sa_k(meta_no_out), N = n_used,
      d = round(meta_no_out$TE.random, 3), ci_lower = round(meta_no_out$lower.random, 3),
      ci_upper = round(meta_no_out$upper.random, 3), I2 = round(meta_no_out$I2 * 100, 1),
      tau2 = round(meta_no_out$tau2, 4), p_value = round(meta_no_out$pval.random, 4),
      n_outliers_excluded = n_excl, status = "completed"
    )
  }
  bind_rows(rows)
}

# =============================================================================
# STAGE 3 — driver: run all 6 analyses, write the joint summaries + run log
# =============================================================================
results_dir <- file.path(PROJECT_ROOT, "results")
region_dir  <- file.path(results_dir, "N_by_region")
country_dir <- file.path(results_dir, "N_by_country")
joint_dir   <- file.path(results_dir, "joint_summary")
dir.create(region_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(country_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(joint_dir, showWarnings = FALSE, recursive = TRUE)

cat(strrep("=", 70), "\n")
cat("STAGE 1: load + validate pool\n")
cat(strrep("=", 70), "\n")
pool_full <- load_and_validate_pool()
save_derived_pool(pool_full)
cat("Loaded and normalized", nrow(pool_full), "rows.\n\n")

DOMAIN_KEYS <- names(DOMAIN_SPEC)

results_by_domain <- list()
for (dk in DOMAIN_KEYS) {
  cat(strrep("=", 70), "\n")
  cat("ANALYSIS:", DOMAIN_SPEC[[dk]]$display, "(", dk, ")\n")
  cat(strrep("=", 70), "\n")
  save_path <- file.path(results_dir, dk)
  results_by_domain[[dk]] <- run_domain_ma(pool_full, dk, save_path, region_dir, country_dir)
  cat("\n")
}

# ── Joint summary across the six analyses ────────────────────────────────────
joint_rows <- lapply(DOMAIN_KEYS, function(dk) {
  m <- results_by_domain[[dk]]
  label <- DOMAIN_SPEC[[dk]]$display
  if (is.null(m)) {
    data.frame(analysis = label, k = NA_integer_, N = NA_integer_, d = NA_real_,
               ci_lower = NA_real_, ci_upper = NA_real_, I2 = NA_real_, tau2 = NA_real_,
               p_value = NA_real_, status = "not_estimable")
  } else {
    data.frame(analysis = label, k = sa_k(m), N = sum(m$data$n, na.rm = TRUE),
               d = round(m$TE.random, 3), ci_lower = round(m$lower.random, 3),
               ci_upper = round(m$upper.random, 3), I2 = round(m$I2 * 100, 1),  # I2 stored as fraction (0-1) by meta::metagen
               tau2 = round(m$tau2, 4), p_value = round(m$pval.random, 4), status = "completed")
  }
})
joint_summary <- bind_rows(joint_rows)
write.csv(joint_summary, file.path(joint_dir, "joint_summary_six_analyses.csv"), row.names = FALSE)
cat("\n=== JOINT SUMMARY (six analyses) ===\n"); print(joint_summary)

# ── No-outliers joint summary — same shape, built from each domain's
# outliers-excluded model (see build_no_outliers_summary() above) ──────────
cat("\n", strrep("=", 70), "\nNO-OUTLIERS JOINT SUMMARY\n", strrep("=", 70), "\n", sep = "")
no_outliers_summary <- build_no_outliers_summary(pool_full)
write.csv(no_outliers_summary, file.path(joint_dir, "joint_summary_six_analyses_no_outliers.csv"), row.names = FALSE)
print(no_outliers_summary)

# ── Run log (every procedure attempted, across all domains) ─────────────────
run_log_df <- bind_rows(RUN_LOG$rows)
write.csv(run_log_df, file.path(results_dir, "run_log_procedures.csv"), row.names = FALSE)

status_tally <- run_log_df %>% count(status) %>% arrange(desc(n))
cat("\n=== RUN LOG STATUS TALLY (", nrow(run_log_df), "procedures attempted) ===\n")
print(status_tally)

errors <- run_log_df %>% filter(status == "error")
if (nrow(errors) > 0) {
  cat("\n=== PROCEDURES THAT ERRORED (", nrow(errors), ") ===\n")
  print(errors %>% select(domain, analysis, moderator, message))
}

# ── Results workbook (joint summary + no-outliers summary + run log) ────────
write_xlsx(list(joint_summary = joint_summary, no_outliers_summary = no_outliers_summary, run_log = run_log_df),
           file.path(joint_dir, "Full-results_six_analyses.xlsx"))

cat("\n", strrep("=", 70), "\n", sep = "")
cat("DONE. Results under:", results_dir, "\n")
cat("Run log:", file.path(results_dir, "run_log_procedures.csv"), "\n")
cat("Joint summary:", file.path(joint_dir, "joint_summary_six_analyses.csv"), "\n")
cat(strrep("=", 70), "\n")
