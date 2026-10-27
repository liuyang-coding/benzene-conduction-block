# Author: Yang Liu
# Date: 2026-10-27 | R 4.5.2
# Packages: survival, CMAverse, ggplot2.
# Input: prepared data with Benzene, All_CB (0/1), follow-up days and covariates.

# 1. Prepared omics data
# p <- nmr_prepare(data, metabolites, measured_metabolites)
# 313 features: 249 measured + 64 ratios; median +/-4 IQR, natural log, then Z.
# Continuous covariates are pre-imputed; categorical missingness is a separate level.
nmr_prepare <- function(data, metabolites, measured_metabolites) {
  if (length(metabolites) != 313L || length(measured_metabolites) != 249L ||
      anyDuplicated(measured_metabolites) || !all(measured_metabolites %in% metabolites))
    stop("Expected 313 features and 249 measured columns.")
  stopifnot(all(c("Benzene", "All_CB", "All_CB_followup_days",
    nmr_covariates(), metabolites) %in% names(data)), !anyDuplicated(names(data)),
    !anyDuplicated(metabolites), all(data$All_CB %in% c(0, 1)),
    all(is.finite(data$All_CB_followup_days) & data$All_CB_followup_days > 0),
    all(vapply(data[c("Benzene", nmr_covariates()[1:7])],
      function(x) is.numeric(x) && all(is.finite(x)), logical(1))))
  for (v in nmr_covariates()[8:18]) {
    x <- as.character(data[[v]])
    x[is.na(x) | x == ""] <- "Missing"
    data[[v]] <- factor(x)
  }
  keep <- rowSums(!is.na(data[measured_metabolites])) > 0L
  data <- data[keep, , drop = FALSE]
  zcols <- paste0(nmr_safe(metabolites), "_z")
  if (anyDuplicated(zcols) || any(zcols %in% names(data))) stop("Standardized names collide; supply untransformed inputs.")
  for (i in seq_along(metabolites)) {
    x <- data[[metabolites[i]]]
    med <- stats::median(x, na.rm = TRUE); spread <- stats::IQR(x, na.rm = TRUE)
    outlier <- if (is.finite(spread) && spread > 0) !is.na(x) & (x < med - 4*spread | x > med + 4*spread) else rep(FALSE, length(x))
    x[outlier] <- NA_real_
    nonpositive <- !is.na(x) & x <= 0
    x[nonpositive] <- NA_real_
    data[[zcols[i]]] <- nmr_z(log(x))
  }
  exposure_sd <- stats::sd(data$Benzene)
  if (!is.finite(exposure_sd) || exposure_sd <= 0) stop("Benzene has no positive SD in the NMR sample.")
  # Standardize within the omics sample; missing omics values remain missing.
  data$Benzene_Z <- nmr_z(data$Benzene)
  list(data = data, features = metabolites, z_columns = zcols, family_size = 313L)
}

# 2. Full-spectrum screening
# s <- nmr_screen(p)
# Separate adjusted linear and Cox paths; Cox screening omits Benzene.
nmr_screen <- function(prepared) {
  data <- prepared$data
  features <- prepared$features
  zcols <- prepared$z_columns
  rows_a <- rows_b <- vector("list", length(features))
  for (i in seq_along(features)) {
    z <- zcols[i]
    a <- data.frame(feature = features[i], mediator = z, n = NA_integer_,
      estimate = NA_real_, std.error = NA_real_, statistic = NA_real_,
      p.value = NA_real_, note = NA_character_)
    b <- data.frame(feature = features[i], mediator = z, analysis_type = "All_CB",
      n = NA_integer_, events = NA_integer_, cox_coef = NA_real_, cox_se = NA_real_,
      cox_z = NA_real_, cox_HR = NA_real_, cox_p = NA_real_, note = NA_character_)
    if (sum(!is.na(data[[z]])) <= 10L) {
      a$note <- b$note <- "Not fitted: <=10 nonmissing standardized values"
    } else {
      a <- tryCatch({
        fit <- stats::lm(stats::reformulate(c("Benzene_Z", nmr_covariates()), response = z),
                         data = data, na.action = stats::na.omit)
        cf <- summary(fit)$coefficients["Benzene_Z", ]
        a$n <- stats::nobs(fit)
        a[1, c("estimate", "std.error", "statistic", "p.value")] <- cf
        a
      }, error = function(e) { a$note <- conditionMessage(e); a })
      b <- tryCatch({
        f <- stats::reformulate(c(z, nmr_covariates()),
          response = "survival::Surv(All_CB_followup_days, All_CB)")
        fit <- survival::coxph(f, data = data, na.action = stats::na.omit)
        cf <- summary(fit)$coefficients[z, ]
        b$n <- fit$n; b$events <- fit$nevent
        b$cox_coef <- cf["coef"]; b$cox_se <- cf["se(coef)"]
        b$cox_z <- cf["z"]; b$cox_HR <- exp(cf["coef"]); b$cox_p <- cf["Pr(>|z|)"]
        b
      }, error = function(e) { b$note <- conditionMessage(e); b })
    }
    rows_a[[i]] <- a; rows_b[[i]] <- b
  }
  path_a <- do.call(rbind, rows_a); path_b <- do.call(rbind, rows_b)
  path_a$p_bonferroni <- stats::p.adjust(path_a$p.value, "bonferroni", n = prepared$family_size)
  path_b$cox_p_bonferroni <- stats::p.adjust(path_b$cox_p, "bonferroni", n = prepared$family_size)
  path_a$beta_low <- path_a$estimate - 1.96 * path_a$std.error
  path_a$beta_high <- path_a$estimate + 1.96 * path_a$std.error
  path_b$HR_low <- exp(path_b$cox_coef - 1.96 * path_b$cox_se)
  path_b$HR_high <- exp(path_b$cox_coef + 1.96 * path_b$cox_se)
  selected_a <- path_a$feature[!is.na(path_a$p_bonferroni) & path_a$p_bonferroni < 0.05]
  selected_b <- path_b$feature[!is.na(path_b$cox_p_bonferroni) & path_b$cox_p_bonferroni < 0.05]
  list(path_a = path_a, path_b = path_b, candidates = sort(intersect(selected_a, selected_b)),
       family_size = prepared$family_size)
}

# 3. Mediation
# m <- nmr_mediate(p, s$candidates)
# CMAverse ratio effects are already exponentiated; PM is an untruncated proportion.
nmr_mediate <- function(prepared, candidates) {
  if (!length(candidates) || anyDuplicated(candidates) || !all(candidates %in% prepared$features))
    stop("Invalid intersection candidates.")
  effects <- c("Rte", "Rpnde", "Rtnie", "pm")
  effect_cols <- unlist(lapply(effects, function(x) paste0(x, c("", "_low", "_high", "_p"))))
  rows <- lapply(candidates, function(feature) {
    z <- prepared$z_columns[match(feature, prepared$features)]
    vars <- c("All_CB_followup_days", "All_CB", "Benzene_Z", z, nmr_covariates())
    d <- prepared$data[stats::complete.cases(prepared$data[vars]), vars, drop = FALSE]
    d <- droplevels(d)
    row <- data.frame(analysis_type = "All_CB", mediator_orig = feature, mediator = z,
      n = nrow(d), events = sum(d$All_CB == 1), stringsAsFactors = FALSE)
    for (v in effect_cols) row[[v]] <- NA_real_
    row$inference_method <- "CMAverse paramfunc with delta method"
    row$note <- NA_character_
    tryCatch({
      fit <- CMAverse::cmest(data = d, model = "rb", outcome = "All_CB_followup_days",
        event = "All_CB", exposure = "Benzene_Z", mediator = z, basec = nmr_covariates(),
        EMint = TRUE, mreg = list("linear"), yreg = "coxph", astar = 0, a = 1,
        mval = list(0), estimation = "paramfunc", inference = "delta")
      values <- nmr_extract_cmest(fit)
      for (v in names(values)) row[[v]] <- values[[v]]
      row
    }, error = function(e) { row$note <- conditionMessage(e); row })
  })
  out <- do.call(rbind, rows)
  # BH denominator includes every candidate, including failed fits.
  out$Rtnie_p_fdr <- stats::p.adjust(out$Rtnie_p, "BH", n = length(candidates))
  if (any(!is.na(out$note))) warning("Mediation failures retained in note; resolve before final summary.")
  out
}

# 4. Indirect-effect BH correction and summary
# summary <- nmr_mediation_summary(m, categories)
nmr_mediation_summary <- function(results, categories = NULL) {
  if (!nrow(results) || anyDuplicated(results$mediator_orig) ||
      !all(results$analysis_type == "All_CB")) stop("Expected one All_CB row per candidate.")
  if (any(!is.na(results$note))) stop("Resolve failed mediation models before producing the summary.")
  results$ACME_p_fdr <- stats::p.adjust(results$Rtnie_p, "BH", n = nrow(results))
  results$significant <- !is.na(results$ACME_p_fdr) & results$ACME_p_fdr < 0.05
  results$category <- if (is.null(categories)) "Unspecified" else unname(categories[results$mediator_orig])
  results$category[is.na(results$category) | results$category == ""] <- "Unmatched"
  by_category <- do.call(rbind, lapply(split(results, results$category), function(d) {
    x <- d$pm[d$significant]
    data.frame(category = d$category[1], n_total = nrow(d), n_sig_fdr = length(x),
      pct_sig = 100 * length(x)/nrow(d), mean_PM = if (length(x)) mean(x) else NA_real_,
      median_PM = if (length(x)) stats::median(x) else NA_real_,
      max_PM = if (length(x)) max(x) else NA_real_)
  }))
  list(full = results, significant = results[results$significant, , drop = FALSE], summary = by_category)
}

# 5. Association plot
# nmr_plot_screen(s, path = "a")
nmr_plot_screen <- function(screen, path = c("a", "b")) {
  path <- match.arg(path)
  d <- if (path == "a") screen$path_a else screen$path_b
  d$effect <- if (path == "a") d$estimate else d$cox_coef
  d$raw_p <- if (path == "a") d$p.value else d$cox_p
  d$selected <- !is.na(d$raw_p) & d$raw_p < 0.05/screen$family_size
  d$neg_log10_p <- -log10(pmax(d$raw_p, .Machine$double.xmin))
  ggplot2::ggplot(d, ggplot2::aes(effect, neg_log10_p, colour = selected)) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = -log10(0.05/screen$family_size), linetype = 2) +
    ggplot2::labs(x = if (path == "a") "Standardized beta" else "Cox log HR per feature SD",
                  y = "-log10(raw P)") + ggplot2::theme_classic()
}

# 6. Spearman correlations
# correlations <- nmr_spearman(p)
# Pairwise-complete transformed values and exact pairwise sample sizes.
nmr_spearman <- function(prepared) {
  x <- as.matrix(prepared$data[prepared$z_columns])
  r <- stats::cor(x, method = "spearman", use = "pairwise.complete.obs")
  n <- crossprod(1L * !is.na(x))
  dimnames(r) <- dimnames(n) <- list(prepared$features, prepared$features)
  list(correlation = r, pairwise_n = n)
}

# Helpers
nmr_covariates <- function() c(
  "age", "BMI", "TDI", "sedentary_time", "passive_smoke_home",
  "outdoor_summer", "outdoor_winter", "sex", "ethnicity_2g", "education_2g",
  "drink_status_2g", "smoke_status_2g", "PA_4g", "diet_three",
  "baseline_hypertension", "baseline_dyslipidemia", "baseline_diabetes",
  "employment_3g"
)

nmr_safe <- function(x) {
  out <- gsub("[^A-Za-z0-9_]", ".", x)
  ifelse(grepl("^[0-9]", out), paste0("m", out), out)
}

nmr_z <- function(x) {
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

nmr_extract_cmest <- function(fit) {
  effects <- as.data.frame(summary(fit)$summarydf)
  keys <- c("Rte", "Rpnde", "Rtnie", "pm")
  cols <- c("Estimate", "95% CIL", "95% CIU")
  if (!all(keys %in% rownames(effects)) || !all(cols %in% names(effects)) ||
      any(!is.finite(as.matrix(effects[keys, cols])))) stop("Invalid CMAverse effect table.")
  pcol <- intersect(c("P.val", "P.value"), names(effects))
  out <- list()
  for (k in keys) {
    out[[k]] <- effects[k, "Estimate"]
    out[[paste0(k, "_low")]] <- effects[k, "95% CIL"]
    out[[paste0(k, "_high")]] <- effects[k, "95% CIU"]
    out[[paste0(k, "_p")]] <- if (length(pcol)) effects[k, pcol[1]] else NA_real_
  }
  out
}
