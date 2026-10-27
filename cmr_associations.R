# CMR associations
# Author: Yang Liu | Date: 2026-10-27 | R: 4.5.2
# Inputs: eligible cohort with eid, Benzene and 18 covariates;
# phenotype table with eid and p24100_i2-p24181_i2. Dependencies: base R, stats.

# Benzene-CMR: common complete sample and 246-test BH family ----
cmr_benzene_associations <- function(cohort, phenotypes) {
  adjustment <- cmr_covariates()
  if ("baseline_CB" %in% names(cohort)) stopifnot(all(cohort$baseline_CB == 0))
  data <- cmr_join(cohort[, c("eid", "Benzene", adjustment)], phenotypes)
  traits <- unname(cmr_field_map())
  data <- data[stats::complete.cases(data[, traits, drop = FALSE]) &
                 is.finite(data$Benzene), , drop = FALSE]
  stopifnot(nrow(data) > 0L, all(vapply(data[traits], function(x) all(is.finite(x)), logical(1))))
  cuts <- stats::quantile(data$Benzene, c(1/3, 2/3), names = FALSE, type = 7)
  stopifnot(stats::sd(data$Benzene) > 0, length(unique(cuts)) == 2L)
  data$Benzene_Z <- as.numeric(scale(data$Benzene))
  data$Benzene_3g <- factor(ifelse(data$Benzene <= cuts[1L], "T1",
    ifelse(data$Benzene <= cuts[2L], "T2", "T3")), levels = c("T1", "T2", "T3"))
  for (column in adjustment[1:7]) {
    median_value <- stats::median(data[[column]], na.rm = TRUE)
    if (!is.finite(median_value)) stop("No finite median for ", column)
    data[[column]][is.na(data[[column]])] <- median_value
  }
  for (column in adjustment[8:18]) {
    values <- as.character(data[[column]])
    values[is.na(values) | values == ""] <- "Missing"
    data[[column]] <- factor(values)
  }
  rows <- lapply(traits, function(trait) {
    grouped <- stats::lm(stats::reformulate(c("Benzene_3g", adjustment), response = trait),
                         data = data, na.action = stats::na.fail)
    continuous <- stats::lm(stats::reformulate(c("Benzene_Z", adjustment), response = trait),
                            data = data, na.action = stats::na.fail)
    values <- rbind(cmr_coefficient(grouped, "Benzene_3gT2"),
                    cmr_coefficient(grouped, "Benzene_3gT3"),
                    cmr_coefficient(continuous, "Benzene_Z"))
    data.frame(trait = trait, contrast = c("T2 versus T1", "T3 versus T1", "Per SD"),
               n = nrow(data), values, row.names = NULL, check.names = FALSE)
  })
  results <- do.call(rbind, rows)
  if (nrow(results) != 246L || any(!is.finite(results$p_value))) {
    stop("All 246 planned benzene-CMR tests must be estimable; inspect sample size/collinearity.")
  }
  results$BH_FDR <- stats::p.adjust(results$p_value, method = "BH", n = 246L)
  # Preserve raw units/signs; optional display columns follow the strain convention.
  reverse <- grepl("^(CS_|LS_)", results$trait)
  results$display_beta <- ifelse(reverse, -results$beta, results$beta)
  results$display_ci_low <- ifelse(reverse, -results$ci_high, results$ci_low)
  results$display_ci_high <- ifelse(reverse, -results$ci_low, results$ci_high)
  results
}
# results <- cmr_benzene_associations(cohort, phenotypes)

# Candidate omics-CMR: upstream Z scores and pairwise complete cases ----
# candidates has feature/column; clinical covariates are already prepared.
cmr_omics_associations <- function(omics, phenotypes, candidates, module = c("NMR", "Olink")) {
  module <- match.arg(module)
  adjustment <- cmr_covariates()
  stopifnot(nrow(candidates) > 0L, !anyNA(candidates[, c("feature", "column")]),
            !anyDuplicated(candidates$feature), !anyDuplicated(candidates$column))
  absent <- !candidates$column %in% names(omics)
  if (any(absent)) {
    if (module == "NMR") stop("NMR candidate columns are absent from omics.")
    warning("Olink candidates absent from omics were omitted: ", paste(candidates$feature[absent], collapse = ", "))
    candidates <- candidates[!absent, , drop = FALSE]
  }
  if (!nrow(candidates)) stop("No candidate columns are available.")
  data <- cmr_join(omics[, c("eid", adjustment, as.character(candidates$column))], phenotypes)
  for (column in adjustment[8:18]) data[[column]] <- factor(data[[column]])
  traits <- unname(cmr_field_map())
  rows <- list()
  for (trait in traits) for (i in seq_len(nrow(candidates))) {
    current <- data[, c(trait, as.character(candidates$column[i]), adjustment), drop = FALSE]
    names(current)[1:2] <- c("outcome", "predictor")
    current <- current[stats::complete.cases(current), , drop = FALSE]
    estimate <- c(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
    status <- "insufficient sample or variation"
    if (nrow(current) >= 30L && length(unique(current$outcome)) >= 2L &&
        (module == "NMR" || length(unique(current$predictor)) >= 2L)) {
      fit <- tryCatch(stats::lm(stats::reformulate(c("predictor", adjustment), response = "outcome"),
                               data = current, na.action = stats::na.fail), error = function(e) e)
      if (inherits(fit, "error")) status <- conditionMessage(fit) else {
        estimate <- cmr_coefficient(fit, "predictor")
        status <- if (is.finite(estimate["p_value"])) "tested" else "nonestimable coefficient"
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(module = module, feature = candidates$feature[i],
      trait = trait, n = nrow(current), as.list(estimate), status = status, check.names = FALSE)
  }
  results <- do.call(rbind, rows)
  # Separate correction families; retain each module's failed-test accounting.
  results$BH_FDR <- if (module == "NMR") {
    stats::p.adjust(results$p_value, method = "BH", n = nrow(results))
  } else stats::p.adjust(results$p_value, method = "BH")
  results
}
# results <- cmr_omics_associations(omics_ready, phenotypes, candidates, "NMR")

# Model setup ----
cmr_covariates <- function() c(
  "age", "BMI", "TDI", "sedentary_time", "passive_smoke_home", "outdoor_summer",
  "outdoor_winter", "sex", "ethnicity_2g", "education_2g", "drink_status_2g",
  "smoke_status_2g", "PA_4g", "diet_three", "baseline_hypertension",
  "baseline_dyslipidemia", "baseline_diabetes", "employment_3g"
)

cmr_field_map <- function() {
  labels <- c("LVEDV", "LVESV", "LVSV", "LVEF", "LVCO", "LVmass", "RVEDV",
              "RVESV", "RVSV", "RVEF", "LAVmax", "LAVmin", "LASV", "LAEF",
              "RAVmax", "RAVmin", "RASV", "RAEF", "AAmax", "AAmin", "AAdist",
              "DAmax", "DAmin", "DAdist", paste0("WT_AHA", 1:16), "WT_global",
              paste0("CS_AHA", 1:16), "CS_global", paste0("RS_AHA", 1:16),
              "RS_global", paste0("LS_Seg", 1:6), "LS_global")
  stats::setNames(labels, paste0("p", 24100:24181, "_i2"))
}

cmr_join <- function(participants, phenotypes) {
  mapping <- cmr_field_map()
  stopifnot(!anyNA(participants$eid), !anyDuplicated(participants$eid),
            !anyNA(phenotypes$eid), !anyDuplicated(phenotypes$eid))
  selected <- phenotypes[, c("eid", names(mapping)), drop = FALSE]
  names(selected) <- c("eid", unname(mapping))
  merge(participants, selected, by = "eid", all = FALSE, sort = FALSE)
}

cmr_coefficient <- function(fit, term) {
  coefficients <- summary(fit)$coefficients
  if (!term %in% rownames(coefficients) || fit$df.residual <= 0L) {
    return(c(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_))
  }
  interval <- stats::confint(fit, parm = term, level = 0.95)
  c(beta = unname(coefficients[term, "Estimate"]),
    se = unname(coefficients[term, "Std. Error"]),
    ci_low = unname(interval[1L]), ci_high = unname(interval[2L]),
    p_value = unname(coefficients[term, "Pr(>|t|)"]))
}
