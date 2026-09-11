# =============================================================================
# META-ANÁLISIS — Cohen's d pool, un engine por dominio de outcome
# =============================================================================
# Corre el mismo pipeline (modelo base, forest, outliers, publication bias,
# subgrupos, meta-regresión, sensitivity analyses) una vez sin filtrar ("All")
# y una vez por cada columna de dominio (brain_outcome, mental_health_outcome,
# behavioral_outcome, cognitive_outcome), escribiendo cada corrida a su propia
# carpeta en results/<domain>/.
#
# Reemplaza a Full-MA.R + Mental_Health-MA.R (movidos a legacy/), que eran
# ~300 líneas casi idénticas entre sí — la única diferencia real era un filtro
# de una línea y la carpeta de salida. Ese tipo de duplicación fue la causa
# original del bug de conteo de "mental health" corregido el 2026-07-17 (dos
# lugares calculando el mismo subconjunto de formas distintas). Un solo engine
# parametrizado por dominio elimina esa clase de bug de raíz: cada dominio se
# calcula en un único lugar, de una única forma.
#
# Un estudio con más de un dominio en "yes" entra a más de un MA de dominio —
# no son subsets mutuamente excluyentes, sino "¿este estudio tiene al menos un
# outcome de este dominio?" (mismo criterio que ya usaba Mental_Health-MA.R).
# =============================================================================
rm(list = ls()); gc()

# run_ma.R lives in scripts/, one level below the project root. [Adapted for
# the MA-Brain&Climate public repo layout: the working pipeline this was
# copied from calls this folder "code/"; only this path and the comment were
# changed to match this repo's "scripts/" folder name — no analysis logic
# was touched.]
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_root <- if (length(script_arg) > 0) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
} else {
  normalizePath(getwd())
}

# ── Paquetes ──────────────────────────────────────────────────────────────────
packages <- c("tidyverse", "readxl", "meta", "metafor", "dmetar",
              "esc", "writexl", "scales", "patchwork")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

setwd(project_root)
source(file.path(project_root, "scripts", "domain_utils.R"))

data_dir     <- file.path(project_root, "data")
results_dir  <- file.path(project_root, "results")
region_dir   <- file.path(results_dir, "N_by_region")
dir.create(region_dir, showWarnings = FALSE, recursive = TRUE)
country_dir  <- file.path(results_dir, "N_by_country")
dir.create(country_dir, showWarnings = FALSE, recursive = TRUE)

# ── Carga y normalización compartida (una sola vez para todos los dominios) ──
pool_d_raw <- read.csv(file.path(data_dir, "pool_continuous_d.csv"), stringsAsFactors = FALSE)

normalize_region <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("East Asia", "South Asia", "Southeast Asia")] <- "Asia"
  x
}
normalize_income <- function(x) {
  x <- tolower(trimws(x))
  case_when(
    grepl("high",         x) & !grepl("upper", x) ~ "High income",
    grepl("upper.middle", x)                       ~ "Upper-middle income",
    grepl("middle",       x) & !grepl("lower|upper", x) ~ "Middle income",
    grepl("moderate",     x)                       ~ "Middle income",
    grepl("lower.middle", x)                       ~ "Lower-middle income",
    TRUE                                            ~ x
  )
}
pool_d_raw$region       <- normalize_region(pool_d_raw$region)
pool_d_raw$income_clean <- normalize_income(pool_d_raw$income_status)
pool_d_raw <- prepare_pool_d(pool_d_raw)  # from domain_utils.R: normalises domain flags, adds outcome_domain + study_label

# ── Helpers de salida (idénticos a los de Full-MA.R, ahora con save_path como parámetro) ──
save_txt <- function(save_path, obj, filename) {
  sink(file.path(save_path, filename))
  print(obj)
  sink()
}

make_forest <- function(save_path, model, filename, xlim, xlab) {
  # 2026-08-01 (user request): figures saved as PDF instead of PNG, kept
  # unfiltered/no-outliers forest plots among the converted ones. Old png()
  # call commented out, not deleted, per explicit instruction.
  # png(file.path(save_path, filename), width = 11, height = 9,
  #     units = "in", res = 300, bg = "white")
  pdf(file.path(save_path, sub("\\.png$", ".pdf", filename)), width = 11, height = 9, bg = "white")
  tryCatch(
    # 2026-08-07 (user request): no user-entered plot title (smlab = ""); no
    # study-level SMD/SE columns on the left of the plot (leftcols = studlab
    # only) — meta::forest()'s default leftcols for a metagen object is
    # c("studlab","TE","seTE"), which is where those unwanted columns came from.
    forest(model, sortvar = model$TE, smlab = "", leftcols = c("studlab"),
           prediction = TRUE,
           zero.pval = TRUE, xlim = xlim, xlab = xlab,
           text.random = "Random effects model", text.predict = "Prediction interval",
           addrows.below.overall = 2, print.tau2 = TRUE, print.I2 = TRUE,
           test.overall = TRUE, rightcols = c("ci", "w.random"),
           rightlabs = c("95% CI", "Weight"), col.square = "#2B7CB5",
           col.diamond = "#D94F3D", fontsize = 10),
    error = function(e) message("forest error: ", e$message))
  dev.off()
}

make_funnel <- function(save_path, model, filename) {
  col.contour <- c("gray75", "gray85", "gray95")
  # png(file.path(save_path, filename), width = 8, height = 7, units = "in", res = 300, bg = "white")
  pdf(file.path(save_path, sub("\\.png$", ".pdf", filename)), width = 8, height = 7, bg = "white")
  funnel(model, yaxs = "i", lwd = 1, cex = 1.2, contour = c(0.9, 0.95, 0.99), col.contour = col.contour)
  legend("topright", legend = c("p < 0.1", "p < 0.05", "p < 0.01"), fill = col.contour, bty = "n")
  dev.off()
}

# GOSH can be slow on large pools; kept as a best-effort diagnostic (tryCatch'd
# at the call site), same behaviour as the original Full-MA.R.
run_gosh <- function(save_path, model, filename, xlim = c(-3, 3)) {
  rma_obj <- rma(yi = model$TE, sei = model$seTE, method = model$method.tau, test = "knha")
  res_gosh <- gosh(rma_obj, parallel = "multicore", ncpus = 2, progbar = FALSE)
  # 2026-08-03 (user request, reverted): GOSH has thousands of points, so as a
  # vector PDF it exploded to ~55MB per domain (vs ~115KB as PNG). Kept as PNG
  # while every other figure in this pipeline moved to PDF — see the skipped
  # pdf() line below, left for reference.
  # pdf(file.path(save_path, sub("\\.png$", ".pdf", filename)), width = 8, height = 7, bg = "white")
  png(file.path(save_path, filename), width = 8, height = 7, units = "in", res = 300, bg = "white")
  plot(res_gosh, alpha = 0.4, het = "I2", col = c("gray70", "#2B7CB5"), xlim = xlim)
  dev.off()
}

run_subgroup <- function(save_path, base_model, subgroup_var, label, filename_prefix, format = "pdf") {
  # 2026-08-01 (user request): all subgroup forest plots move to PDF EXCEPT
  # "design" (Observational/RCT/etc.), which stays PNG per explicit
  # instruction — called with format = "png" at that one call site below.
  model_sg <- tryCatch(update(base_model, subgroup = subgroup_var),
    error = function(e) { message("subgroup error (", label, "): ", e$message); NULL })
  if (is.null(model_sg)) return(invisible(NULL))
  save_txt(save_path, model_sg, paste0(filename_prefix, "_subgroup_", label, ".txt"))
  plot_file <- file.path(save_path, paste0(filename_prefix, "_subgroup_", label, ".", format))
  plot_height <- max(8, 0.3 * nrow(base_model$data) + 4)
  if (format == "png") {
    png(plot_file, width = 13, height = plot_height, units = "in", res = 300, bg = "white")
  } else {
    pdf(plot_file, width = 13, height = plot_height, bg = "white")
  }
  tryCatch(forest(model_sg, sortvar = model_sg$TE, prediction = TRUE,
    print.tau2 = TRUE, print.I2 = TRUE, test.subgroup = TRUE,
    leftcols = c("studlab"), rightcols = c("ci", "w.random"), rightlabs = c("95% CI", "Weight"),
    col.square = "#2B7CB5", col.diamond = "#D94F3D", fontsize = 9),
    error = function(e) message("subgroup forest error: ", e$message))
  dev.off()
  invisible(model_sg)
}

run_metareg <- function(save_path, data, formula_str, label, filename_prefix) {
  rma_obj <- tryCatch(
    rma(yi = data$d, sei = data$se_d, mods = as.formula(formula_str),
        method = "REML", test = "knha", data = data),
    error = function(e) { message("Meta-reg error (", label, "): ", e$message); NULL })
  if (is.null(rma_obj)) return(invisible(NULL))
  sink(file.path(save_path, paste0(filename_prefix, "_metareg_", label, ".txt")))
  cat("Meta-regresión:", label, "\nFórmula:", formula_str, "\n\n")
  print(summary(rma_obj))
  cat("\nR²:", round(rma_obj$R2, 2), "%\n")
  sink()
  if (label == "log_n") {
    tryCatch({
      # png(file.path(save_path, paste0(filename_prefix, "_bubble_log_n.png")),
      #     width = 7, height = 6, units = "in", res = 300, bg = "white")
      pdf(file.path(save_path, paste0(filename_prefix, "_bubble_log_n.pdf")),
          width = 7, height = 6, bg = "white")
      regplot(rma_obj, mod = 2, xlab = "log(N)", ylab = "Effect size", col = "#2B7CB5", bg = "#A8C8E8", pch = 21)
      dev.off()
    }, error = function(e) { dev.off(); message("regplot error: ", e$message) })
  }
  invisible(rma_obj)
}

sa_val <- function(obj, field) { if (is.null(obj)) return(NA_real_); val <- obj[[field]]; if (is.null(val) || length(val) == 0) NA_real_ else val }
sa_k   <- function(obj) { if (is.null(obj)) return(NA_integer_); k <- obj[["k"]]; if (is.null(k) || length(k) == 0) NA_integer_ else k }

# =============================================================================
# ENGINE: corre el pipeline completo para un dominio (o "All" si domain_col es NULL)
# =============================================================================
run_domain_ma <- function(pool_d_full, domain_col, domain_label) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("DOMAIN:", domain_label, "\n")
  cat(strrep("=", 60), "\n", sep = "")

  save_path <- file.path(results_dir, domain_label)
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

  pool_d <- if (is.null(domain_col)) pool_d_full else pool_d_full %>% filter(.data[[domain_col]] == "yes")
  cat("Outcome rows in this domain:", nrow(pool_d), "\n")

  # ── k papers y N poblacional por estudio (punto 5) ──────────────────────────
  study_summary <- pool_d %>%
    group_by(study_label, region, income_clean) %>%
    summarise(
      k_outcomes        = n(),
      # max(), not first(): see domain_utils.R::aggregate_d_pool for why.
      n_population      = max(N, na.rm = TRUE),
      design_family     = best_design(design_family),
      outcomes_combined = paste(outcome_name, collapse = " // "),
      d_composite       = mean(d_computed, na.rm = TRUE),
      se_composite      = composite_se(d_se),
      .groups = "drop")
  write.csv(study_summary, file.path(save_path, "00_study_summary.csv"), row.names = FALSE)

  pool_d_agg <- aggregate_d_pool(pool_d)
  k <- nrow(pool_d_agg)
  cat("Independent study composites (k):", k, "\n")

  if (k < 2) {
    msg <- sprintf(paste(
      "Domain '%s': only %d independent study composite(s) in the Cohen's d pool",
      "after aggregation — not enough to run a random-effects meta-analysis (needs >= 2).",
      "See 00_study_summary.csv for what is available for this domain.", sep = "\n"),
      domain_label, k)
    writeLines(msg, file.path(save_path, "00_NOTE_insufficient_data.txt"))
    cat(msg, "\n")
    return(invisible(NULL))
  }

  n_by_region <- pool_d_agg %>% group_by(region) %>%
    summarise(k = n(), N = sum(n, na.rm = TRUE), .groups = "drop")
  write.csv(n_by_region, file.path(region_dir, paste0("N_by_region_", domain_label, ".csv")), row.names = FALSE)

  n_by_country <- pool_d_agg %>% group_by(country) %>%
    summarise(k = n(), N = sum(n, na.rm = TRUE), .groups = "drop")
  write.csv(n_by_country, file.path(country_dir, paste0("N_by_country_", domain_label, ".csv")), row.names = FALSE)

  # ── 2. Modelo base — random effects (REML + Hartung-Knapp) ──────────────────
  meta_d <- metagen(
    TE = d, seTE = se_d, studlab = study_label, data = pool_d_agg,
    sm = "SMD", method.tau = "REML", prediction = TRUE,
    common = FALSE, random = TRUE, method.random.ci = "HK",
    title = paste("Meta-análisis pool continuo (SMD) —", domain_label)
  )
  save_txt(save_path, meta_d, "01_meta_d_base.txt")
  make_forest(save_path, meta_d, "02_forest_d.png", c(-3, 3), "Cohen's d")

  # ── 3. Outliers ───────────────────────────────────────────────────────────
  out_d <- tryCatch(find.outliers(meta_d), error = function(e) { message("outliers: ", e$message); NULL })
  meta_d_no_out <- meta_d
  if (!is.null(out_d)) {
    save_txt(save_path, out_d, "03_outliers_d.txt")
    meta_d_no_out <- tryCatch(update(meta_d, exclude = meta_d$studlab %in% out_d$out.study.random),
                               error = function(e) meta_d)
    save_txt(save_path, meta_d_no_out, "03_meta_d_no_outliers.txt")
    make_forest(save_path, meta_d_no_out, "03_forest_d_no_outliers.png", c(-3, 3), "Cohen's d")
  }

  # ── 4. Publication bias ──────────────────────────────────────────────────
  make_funnel(save_path, meta_d, "04_funnel_d.png")
  save_txt(save_path, tryCatch(eggers.test(meta_d), error = function(e) e), "04_egger_d.txt")
  if (!identical(meta_d_no_out, meta_d)) {
    save_txt(save_path, tryCatch(eggers.test(meta_d_no_out), error = function(e) e), "04_egger_d_no_out.txt")
  }
  tryCatch({
    # png(file.path(save_path, "04_pcurve_d.png"), width = 7, height = 7, units = "in", res = 300, bg = "white")
    pdf(file.path(save_path, "04_pcurve_d.pdf"), width = 7, height = 7, bg = "white")
    sink(file.path(save_path, "04_pcurve_d.txt"))
    print(pcurve(meta_d))
    sink()
    dev.off()
  }, error = function(e) { message("P-curve error: ", conditionMessage(e)); try(sink()); try(dev.off()) })

  # ── 5. GOSH ───────────────────────────────────────────────────────────────
  tryCatch(run_gosh(save_path, meta_d, "05_gosh_d.png", xlim = c(-2, 2)),
           error = function(e) message("GOSH error: ", e$message))

  # ── 6. Subgrupos ──────────────────────────────────────────────────────────
  # "design" stays PNG per explicit user request (2026-08-01) — every other
  # figure in this pipeline moved to PDF, this one and the sensitivity
  # comparison forest (08_sa_comparison_forest.png, below) did not.
  run_subgroup(save_path, meta_d, pool_d_agg$design_family, "design", "06_d", format = "png")
  # outcome_domain solo tiene sentido cruzando dominios en el pool sin filtrar;
  # dentro de un dominio ya filtrado, casi todas las filas comparten el mismo
  # outcome_domain por construcción.
  #
  # Recoding (2026-07-31, user decision): "mental_health/behavioral" (rows
  # double-tagged on both flags simultaneously) is folded into "mental_health"
  # — same domain family, not a distinct group. "unclassified" (no domain flag
  # set at all) is dropped from this specific subgroup/metareg only — it isn't
  # a real outcome domain, so comparing effect sizes against a "no domain"
  # bucket is not meaningful. Both were k=3/k=2, too sparse to interpret on
  # their own. This does not change the "All" pool or its k=36 base model —
  # only this domain-subgroup view uses a filtered/recoded copy.
  if (identical(domain_label, "All")) {
    domain_recoded <- dplyr::if_else(pool_d_agg$outcome_domain == "mental_health/behavioral",
                                      "mental_health", pool_d_agg$outcome_domain)
    domain_keep <- domain_recoded != "unclassified"
    pool_d_agg_domain_sg <- pool_d_agg[domain_keep, ]
    meta_d_domain <- metagen(
      TE = d, seTE = se_d, studlab = study_label, data = pool_d_agg_domain_sg,
      sm = "SMD", method.tau = "REML", prediction = TRUE,
      common = FALSE, random = TRUE, method.random.ci = "HK",
      title = paste("Meta-análisis pool continuo (SMD) —", domain_label, "(outcome-domain recoded)")
    )
    run_subgroup(save_path, meta_d_domain, domain_recoded[domain_keep], "outcome_domain", "06_d")
  }
  run_subgroup(save_path, meta_d, pool_d_agg$region, "region", "06_d")
  run_subgroup(save_path, meta_d, pool_d_agg$income_clean, "income", "06_d")

  # ── 7. Meta-regresión ─────────────────────────────────────────────────────
  design_levels <- levels(factor(pool_d_agg$design_family))
  design_ref <- if ("observational" %in% design_levels) "observational" else if ("quasi-experimental" %in% design_levels) "quasi-experimental" else design_levels[1]
  pool_d_agg$design_fct <- relevel(factor(pool_d_agg$design_family), ref = design_ref)

  region_levels <- levels(factor(pool_d_agg$region))
  region_ref <- if ("Europe" %in% region_levels) "Europe" else region_levels[1]
  pool_d_agg$region_fct <- relevel(factor(pool_d_agg$region), ref = region_ref)

  income_levels <- levels(factor(pool_d_agg$income_clean))
  income_ref <- if ("High income" %in% income_levels) "High income" else income_levels[1]
  pool_d_agg$income_fct <- relevel(factor(pool_d_agg$income_clean), ref = income_ref)

  run_metareg(save_path, pool_d_agg, "~ design_fct", "design", "07_d")
  if (identical(domain_label, "All")) {
    # Same recoding as the outcome_domain subgroup above (mental_health/behavioral
    # folded into mental_health; unclassified excluded) — see comment there.
    domain_recoded <- dplyr::if_else(pool_d_agg$outcome_domain == "mental_health/behavioral",
                                      "mental_health", pool_d_agg$outcome_domain)
    domain_keep <- domain_recoded != "unclassified"
    pool_d_agg_domain <- pool_d_agg[domain_keep, ]
    pool_d_agg_domain$outcome_domain_recoded <- domain_recoded[domain_keep]
    domain_levels <- levels(factor(pool_d_agg_domain$outcome_domain_recoded))
    domain_ref <- if ("mental_health" %in% domain_levels) "mental_health" else domain_levels[1]
    pool_d_agg_domain$domain_fct <- relevel(factor(pool_d_agg_domain$outcome_domain_recoded), ref = domain_ref)
    run_metareg(save_path, pool_d_agg_domain, "~ domain_fct", "outcome_domain", "07_d")
  }
  run_metareg(save_path, pool_d_agg, "~ region_fct", "region", "07_d")
  run_metareg(save_path, pool_d_agg, "~ income_fct", "income", "07_d")
  run_metareg(save_path, pool_d_agg, "~ log(n)", "log_n", "07_d")

  # ── 8. Sensitivity analyses ──────────────────────────────────────────────
  # Nota: la SA "mental health only" que Full-MA.R corría aparte (re-agregada
  # a mano para no repetir el bug de conteo del 2026-07-17) ya no hace falta:
  # results/mental_health/ es ahora un MA completo por derecho propio, calculado
  # con esta misma función. Mantener esa SA acá habría vuelto a duplicar la
  # misma lógica en dos lugares — exactamente lo que causó el bug original.
  loo_d <- tryCatch(metainf(meta_d), error = function(e) { message("loo_d: ", e$message); NULL })
  if (!is.null(loo_d)) {
    # png(file.path(save_path, "08_loo_d.png"), width = 10, height = 8, units = "in", res = 300, bg = "white")
    pdf(file.path(save_path, "08_loo_d.pdf"), width = 10, height = 8, bg = "white")
    tryCatch(forest(loo_d, xlab = "Cohen's d", leftcols = c("studlab"), col.square = "#2B7CB5", fontsize = 9),
             error = function(e) message("loo_d forest: ", e$message))
    dev.off()
  }

  meta_d_rct <- tryCatch(update(meta_d, subset = pool_d_agg$design_family %in% c("RCT", "RCT_cluster")),
                          error = function(e) { message("SA RCT: ", e$message); NULL })
  if (!is.null(meta_d_rct)) save_txt(save_path, meta_d_rct, "08_sa_d_rct_only.txt")

  meta_d_obs <- tryCatch(update(meta_d, subset = pool_d_agg$design_family == "observational"),
                          error = function(e) { message("SA obs: ", e$message); NULL })
  if (!is.null(meta_d_obs)) save_txt(save_path, meta_d_obs, "08_sa_d_observational_only.txt")

  studies_d_direct <- pool_d %>% filter(es_pool == "d_direct") %>% pull(study_label) %>% unique()
  meta_d_direct_only <- tryCatch(update(meta_d, subset = pool_d_agg$study_label %in% studies_d_direct),
                                  error = function(e) { message("SA d-direct: ", e$message); NULL })
  if (!is.null(meta_d_direct_only)) save_txt(save_path, meta_d_direct_only, "08_sa_d_direct_only.txt")

  tau_comparison <- data.frame(estimator = c("PM", "REML", "DL", "HS"), tau2_d = NA_real_, I2_d = NA_real_, TE_d = NA_real_)
  for (i in seq_along(tau_comparison$estimator)) {
    m <- tryCatch(update(meta_d, method.tau = tau_comparison$estimator[i], method.random.ci = "classic"), error = function(e) NULL)
    if (!is.null(m)) { tau_comparison$tau2_d[i] <- m$tau2; tau_comparison$I2_d[i] <- m$I2; tau_comparison$TE_d[i] <- m$TE.random }
  }
  write.csv(tau_comparison, file.path(save_path, "08_sa_tau_estimators.csv"), row.names = FALSE)

  sa_summary <- data.frame(
    analysis = c("Base (All studies)", "Without outliers", "RCT only", "Observational only", "Direct Cohen's d only"),
    k        = c(sa_k(meta_d), sa_k(meta_d_no_out), sa_k(meta_d_rct), sa_k(meta_d_obs), sa_k(meta_d_direct_only)),
    d        = round(c(sa_val(meta_d,"TE.random"), sa_val(meta_d_no_out,"TE.random"), sa_val(meta_d_rct,"TE.random"), sa_val(meta_d_obs,"TE.random"), sa_val(meta_d_direct_only,"TE.random")), 3),
    ci_lower = round(c(sa_val(meta_d,"lower.random"), sa_val(meta_d_no_out,"lower.random"), sa_val(meta_d_rct,"lower.random"), sa_val(meta_d_obs,"lower.random"), sa_val(meta_d_direct_only,"lower.random")), 3),
    ci_upper = round(c(sa_val(meta_d,"upper.random"), sa_val(meta_d_no_out,"upper.random"), sa_val(meta_d_rct,"upper.random"), sa_val(meta_d_obs,"upper.random"), sa_val(meta_d_direct_only,"upper.random")), 3),
    I2       = round(c(sa_val(meta_d,"I2"), sa_val(meta_d_no_out,"I2"), sa_val(meta_d_rct,"I2"), sa_val(meta_d_obs,"I2"), sa_val(meta_d_direct_only,"I2")), 1),
    p_value  = round(c(sa_val(meta_d,"pval.random"), sa_val(meta_d_no_out,"pval.random"), sa_val(meta_d_rct,"pval.random"), sa_val(meta_d_obs,"pval.random"), sa_val(meta_d_direct_only,"pval.random")), 4)
  )
  write.csv(sa_summary, file.path(save_path, "08_sensitivity_summary.csv"), row.names = FALSE)
  cat("\n=== SENSITIVITY SUMMARY —", domain_label, "===\n")
  print(sa_summary)

  # Stays PNG per explicit user request (2026-08-01) — this is the sensitivity
  # comparison plot that includes the RCT-only/Observational-only rows, kept
  # out of the PNG->PDF conversion along with 06_d_subgroup_design.png.
  png(file.path(save_path, "08_sa_comparison_forest.png"), width = 9, height = 6, units = "in", res = 300, bg = "white")
  par(mar = c(5, max(nchar(sa_summary$analysis)) * 0.6 + 1, 4, 1))
  plot(x = sa_summary$d, y = nrow(sa_summary):1,
       xlim = c(min(sa_summary$ci_lower, na.rm = TRUE) - 0.1, max(sa_summary$ci_upper, na.rm = TRUE) + 0.1),
       yaxt = "n", xlab = "Cohen's d (RE)", ylab = "", pch = 18, cex = 1.5, col = "#2B7CB5",
       main = "")
  segments(x0 = sa_summary$ci_lower, x1 = sa_summary$ci_upper, y0 = nrow(sa_summary):1, y1 = nrow(sa_summary):1, col = "#2B7CB5", lwd = 2)
  axis(2, at = nrow(sa_summary):1, labels = sa_summary$analysis, las = 1, cex.axis = 0.8)
  abline(v = 0, lty = 2, col = "gray50")
  dev.off()

  cat(sprintf("\n%s\nFINAL SUMMARY — %s\n%s\n\n", strrep("=", 60), domain_label, strrep("=", 60)))
  cat(sprintf("  k = %d studies\n", meta_d$k))
  cat(sprintf("  d = %.3f [%.3f, %.3f]\n", meta_d$TE.random, meta_d$lower.random, meta_d$upper.random))
  cat(sprintf("  I² = %.1f%%   τ² = %.4f\n", meta_d$I2, meta_d$tau2))
  cat(sprintf("  p = %.4f\n", meta_d$pval.random))
  cat(sprintf("  Results saved: %s\n", save_path))

  invisible(meta_d)
}

# =============================================================================
# RUN — "All" (sin subset) + un MA por columna de dominio de outcome
# =============================================================================
DOMAINS <- list(
  All           = NULL,
  mental_health = "mental_health_outcome",
  behavioral    = "behavioral_outcome"
)

results_by_domain <- list()
for (domain_label in names(DOMAINS)) {
  results_by_domain[[domain_label]] <- run_domain_ma(pool_d_raw, DOMAINS[[domain_label]], domain_label)
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("ALL DOMAINS DONE:", paste(names(DOMAINS), collapse = ", "), "\n")
cat(strrep("=", 60), "\n", sep = "")
