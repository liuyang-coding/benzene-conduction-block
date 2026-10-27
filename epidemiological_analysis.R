# Author: Yang Liu
# Date: 2026-10-27
# R version: 4.5.2

# Packages ----
# survival, tableone, rms

# Prepared data ----
# data <- readRDS("data/inputs/prepared_cohort.rds")
# Prepared incident cohort; Benzene is the exposure variable.

# Covariate preparation ----
# model_data <- epi_prepare_covariates(data)
epi_prepare_covariates <- function(data) {
  epi_impute(data)
}

# Baseline characteristics ----
# table1 <- epi_descriptive(data)
epi_descriptive <- function(data, include_prs = TRUE) {
  variables <- epi_variables()
  for (v in variables$categorical) {
    data[[v]] <- factor(data[[v]])
  }
  if (include_prs) {
    stopifnot(all(is.finite(data$PRS)))
    data$Genetic_risk_category <- epi_tertiles(data$PRS)$group
    variables$all <- c(variables$all, "Genetic_risk_category")
    variables$categorical <- c(variables$categorical, "Genetic_risk_category")
  }
  tableone::CreateTableOne(vars = variables$all, strata = "All_CB", data = data,
    factorVars = variables$categorical, includeNA = TRUE, test = TRUE,
    smd = FALSE, addOverall = TRUE)
}

# Cox regression ----
# cox <- epi_cox(data)
epi_cox <- function(data, outcomes = c("All_CB", "AB", "LB", "RB", "Other_CB")) {
  adjusted15 <- c("age", "sex", "TDI", "ethnicity_2g", "education_2g", "employment_3g",
    "BMI", "drink_status_2g", "smoke_status_2g", "sedentary_time", "PA_4g",
    "outdoor_summer", "outdoor_winter", "diet_three", "passive_smoke_home")
  adjustments <- list("Model 1" = character(), "Model 2" = adjusted15,
    "Fully adjusted" = c(adjusted15, "baseline_hypertension", "baseline_dyslipidemia",
                         "baseline_diabetes"))
  fits <- list(); rows <- list(); diagnostics <- list(); audits <- list()
  for (outcome in outcomes) {
    d <- epi_prepare(data, outcome)
    cuts <- attr(d, "tertile_cuts"); sd_x <- attr(d, "exposure_sd")
    audits[[outcome]] <- cbind(outcome = outcome, epi_counts(d), SD = sd_x,
      IQR = stats::IQR(d$exposure), T1_upper = cuts[1], T2_upper = cuts[2])
    for (model_name in names(adjustments)) {
      group_fits <- lapply(c("exposure_SD", "exposure_tertile", "exposure_score"),
        function(term) epi_fit(d, term, adjustments[[model_name]]))
      names(group_fits) <- c("continuous", "tertile", "trend")
      trend <- epi_extract(group_fits$trend, "exposure_score")
      results <- rbind(
        cbind(contrast = "Per SD increase", epi_extract(group_fits$continuous, "exposure_SD")),
        cbind(contrast = c("T2 vs T1", "T3 vs T1"),
              epi_extract(group_fits$tertile, c("exposure_tertileT2", "exposure_tertileT3"))))
      results$P_trend <- c(NA_real_, rep(trend$P, 2))
      results$log_P_trend <- c(NA_real_, rep(trend$log_P, 2))
      key <- paste(outcome, model_name, sep = "/")
      rows[[key]] <- cbind(exposure = "Benzene", outcome = outcome, model = model_name,
        results, participants = nrow(d), cases = sum(d$event), SD = sd_x)
      fits[[key]] <- group_fits
      diagnostics[[key]] <- lapply(group_fits[c("continuous", "tertile")], survival::cox.zph)
    }
  }
  results <- do.call(rbind, rows); rownames(results) <- NULL
  results$P_BH_within_exposure <- ifelse(results$outcome == "All_CB", NA_real_, results$P)
  list(results = results, sample = do.call(rbind, audits), models = fits, PH = diagnostics)
}

# E-value input ----
# hr_table <- epi_evalue_input(cox)
epi_evalue_input <- function(cox_result) {
  d <- cox_result$results
  d <- d[d$model=="Fully adjusted",c("exposure","outcome","contrast","HR","lower","upper",
                                     "participants","cases"),drop=FALSE]
  names(d)[match(c("HR","lower","upper"),names(d))] <- c("hr","lcl","ucl")
  stopifnot(nrow(d)==15L,!anyDuplicated(d[c("exposure","outcome","contrast")]))
  rownames(d) <- NULL
  d
}

# Restricted cubic splines ----
# rcs <- epi_rcs(data, outcome = "All_CB")
epi_rcs <- function(data, outcome = "All_CB") {
  d <- epi_impute(data); d$x <- d$Benzene
  limits <- stats::quantile(d$x, c(.025, .975), na.rm = TRUE, type = 7)
  d <- d[is.finite(d$x) & d$x >= limits[1] & d$x <= limits[2], , drop = FALSE]
  reference <- stats::median(d$x)
  schemes <- list(`3` = c(.10,.50,.90), `4` = c(.05,.35,.65,.95),
                  `5` = c(.05,.275,.50,.725,.95))
  knots <- lapply(schemes, function(p) unname(stats::quantile(d$x, p, type = 7)))
  d$time <- d[[paste0(outcome, "_followup_days")]] / 365.25; d$event <- d[[outcome]]
  d <- droplevels(d[is.finite(d$time) & d$time > 0 & !is.na(d$event), , drop = FALSE])
  stopifnot(all(d$event %in% 0:1), sum(d$event) > 0)
  rcs <- rms::rcs
  fits <- lapply(knots, function(k) {
    stopifnot(length(unique(k)) == length(k))
    f <- epi_formula(paste0("rcs(x, c(", paste(sprintf("%.17g", k), collapse=","), "))"))
    environment(f) <- environment()
    fit <- rms::cph(f, data = d, x = TRUE, y = TRUE, surv = TRUE,
                    method = "efron", na.action = stats::na.fail)
    stopifnot(!isTRUE(fit$fail), all(is.finite(stats::coef(fit))))
    fit
  })
  bic <- do.call(rbind, lapply(names(fits), function(k) {
    fit <- fits[[k]]; ll <- tail(fit$loglik, 1); p <- length(stats::coef(fit))
    data.frame(knots = as.integer(k), log_partial_likelihood = ll, parameters = p,
      events = sum(d$event), BIC = -2*ll + p*log(sum(d$event)))
  }))
  fit <- fits[["3"]]
  grid <- sort(unique(c(seq(min(d$x), max(d$x), length.out=200), reference)))
  newdata <- d[rep(1, length(grid)), , drop = FALSE]; newdata$x <- grid
  design <- stats::predict(fit, newdata, type = "x")
  contrast <- sweep(design, 2, design[match(reference, grid), ], "-")
  beta <- drop(contrast %*% stats::coef(fit))
  se <- sqrt(pmax(0, rowSums((contrast %*% stats::vcov(fit))*contrast)))
  curve <- data.frame(concentration = grid, HR = exp(beta),
    lower = exp(beta-1.96*se), upper = exp(beta+1.96*se))
  list(curve = curve, anova = stats::anova(fit), PH = survival::cox.zph(fit),
    BIC = bic, knots = knots, reference = reference, trim_limits = limits,
    sample = epi_counts(d), models = fits)
}

# Subgroup interactions ----
# subgroups <- epi_subgroup(data)
epi_subgroup <- function(data) {
  medians <- vapply(data[epi_variables()$continuous], stats::median, numeric(1), na.rm=TRUE)
  d <- epi_impute(data)
  d$age_group <- factor(as.integer(d$age >= 60), levels=0:1)
  d$BMI_2g <- factor(as.integer(d$BMI >= 30), levels=0:1)
  for (v in c("TDI", "sedentary_time", "passive_smoke_home", "outdoor_summer", "outdoor_winter")) {
    high <- if (v %in% c("TDI", "passive_smoke_home")) d[[v]] > medians[v] else d[[v]] >= medians[v]
    d[[paste0(v, "_2g")]] <- factor(as.integer(high), levels=0:1)
  }
  d <- epi_prepare(d, impute=FALSE)
  subgroups <- c("age_group", "sex", "TDI_2g", "ethnicity_2g", "education_2g",
    "employment_3g", "BMI_2g", "drink_status_2g", "smoke_status_2g", "sedentary_time_2g",
    "PA_4g", "outdoor_summer_2g", "outdoor_winter_2g", "diet_three", "passive_smoke_home_2g",
    "baseline_hypertension", "baseline_dyslipidemia", "baseline_diabetes")
  rows <- list(); lrt <- list(); models <- list()
  for (v in subgroups) {
    covariates <- setdiff(epi_variables()$all, v)
    additive <- epi_fit(d, paste("exposure_tertile +", v), covariates)
    interaction <- epi_fit(d, paste("exposure_tertile *", v), covariates)
    statistic <- 2*(as.numeric(stats::logLik(interaction))-as.numeric(stats::logLik(additive)))
    df <- attr(stats::logLik(interaction), "df")-attr(stats::logLik(additive), "df")
    stopifnot(statistic > -1e-6, df > 0)
    lrt[[v]] <- data.frame(subgroup=v, chi_square=max(0,statistic), df=df,
      P=stats::pchisq(max(0,statistic), df, lower.tail=FALSE))
    for (level in levels(d[[v]])) for (group in c("T2", "T3")) {
      group_rows <- as.character(d[[v]]) == level
      newdata <- d[rep(which(group_rows)[1], 2), , drop=FALSE]
      newdata$exposure_tertile <- factor(c("T1", group), levels=c("T1","T2","T3"))
      design <- stats::model.matrix(stats::delete.response(stats::terms(interaction)), newdata,
        contrasts.arg=interaction$contrasts, xlev=interaction$xlevels)
      terms <- names(stats::coef(interaction)); contrast <- design[2,terms]-design[1,terms]
      beta <- sum(contrast*stats::coef(interaction))
      se <- sqrt(drop(t(contrast)%*%stats::vcov(interaction)%*%contrast))
      selected <- group_rows & d$exposure_tertile == group
      rows[[paste(v,level,group)]] <- cbind(subgroup=v, level=level,
        contrast=paste(group,"vs T1"), epi_effect(beta,se), epi_counts(d[selected,,drop=FALSE]))
    }
    models[[v]] <- list(additive=additive, interaction=interaction)
  }
  list(results=do.call(rbind,rows), interaction=do.call(rbind,lrt), models=models,
       tertile_cuts=attr(d,"tertile_cuts"), subgroup_medians=medians)
}

# Helper functions ----
epi_variables <- function() {
  continuous <- c("age", "BMI", "TDI", "sedentary_time", "passive_smoke_home",
                  "outdoor_summer", "outdoor_winter")
  categorical <- c("sex", "ethnicity_2g", "education_2g", "drink_status_2g",
    "smoke_status_2g", "PA_4g", "diet_three", "baseline_hypertension",
    "baseline_dyslipidemia", "baseline_diabetes", "employment_3g")
  list(continuous = continuous, categorical = categorical,
       all = c(continuous, categorical))
}

epi_impute <- function(data) {
  variables <- epi_variables()
  for (v in variables$continuous) {
    med <- stats::median(data[[v]], na.rm = TRUE)
    if (!is.finite(med)) stop("No finite median for ", v)
    data[[v]][is.na(data[[v]])] <- med
  }
  for (v in variables$categorical) {
    z <- as.character(data[[v]])
    z[is.na(z) | z == ""] <- "Missing"
    data[[v]] <- factor(z)
  }
  data
}

epi_tertiles <- function(x, cuts = NULL) {
  if (is.null(cuts)) cuts <- unname(stats::quantile(x, c(1/3, 2/3), type = 7))
  stopifnot(length(cuts) == 2L, all(is.finite(cuts)), cuts[1] < cuts[2])
  list(group = cut(x, c(-Inf, cuts, Inf), labels = c("T1", "T2", "T3"),
                   right = TRUE), cuts = cuts)
}

epi_prepare <- function(data, outcome = "All_CB", impute = TRUE) {
  if (impute) data <- epi_impute(data)
  data$exposure <- data$Benzene
  data$time <- data[[paste0(outcome, "_followup_days")]] / 365.25
  data$event <- data[[outcome]]
  stopifnot(is.numeric(data$event), all(is.na(data$event) | data$event %in% 0:1))
  data <- droplevels(data[is.finite(data$exposure) & is.finite(data$time) &
    data$time > 0 & !is.na(data$event), , drop = FALSE])
  stopifnot(nrow(data) > 1L, sum(data$event) > 0)
  tertiles <- epi_tertiles(data$exposure)
  exposure_sd <- stats::sd(data$exposure)
  stopifnot(is.finite(exposure_sd), exposure_sd > 0)
  data$exposure_tertile <- tertiles$group
  data$exposure_score <- as.integer(data$exposure_tertile)
  data$exposure_SD <- data$exposure / exposure_sd
  attr(data, "tertile_cuts") <- tertiles$cuts
  attr(data, "exposure_sd") <- exposure_sd
  data
}

epi_formula <- function(term, covariates = epi_variables()$all,
                           response = "survival::Surv(time, event)") {
  stats::as.formula(paste(response, "~", paste(c(term, covariates), collapse = " + ")))
}

epi_effect <- function(beta, se, df = Inf) {
  critical <- if (is.finite(df)) stats::qt(.975, df) else stats::qnorm(.975)
  log_p <- if (is.finite(df)) log(2) + stats::pt(abs(beta/se), df,
    lower.tail = FALSE, log.p = TRUE) else stats::pchisq((beta/se)^2, 1,
    lower.tail = FALSE, log.p = TRUE)
  data.frame(estimate = unname(beta), SE = unname(se), df = df, HR = exp(unname(beta)),
    lower = exp(unname(beta) - critical * unname(se)),
    upper = exp(unname(beta) + critical * unname(se)), P = exp(log_p), log_P = log_p)
}

epi_extract <- function(model, terms = names(stats::coef(model))) {
  beta <- stats::coef(model)[terms]
  se <- sqrt(diag(stats::vcov(model)))[terms]
  cbind(term = terms, epi_effect(beta, se), row.names = NULL)
}

epi_fit <- function(data, term, covariates = epi_variables()$all,
                      response = "survival::Surv(time, event)") {
  model <- survival::coxph(epi_formula(term, covariates, response), data = data,
    ties = "efron", x = TRUE, y = TRUE, na.action = stats::na.fail)
  stopifnot(model$n == nrow(data), all(is.finite(stats::coef(model))))
  model
}

epi_counts <- function(data) {
  data.frame(participants = nrow(data), cases = sum(data$event),
    person_years = sum(data$time), median_followup_years = stats::median(data$time),
    incidence_10000_py = sum(data$event)/sum(data$time)*10000)
}
