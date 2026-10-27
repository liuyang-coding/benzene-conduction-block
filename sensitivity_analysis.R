# Author: Yang Liu
# Date: 2026-10-27
# R version: 4.5.2

# Packages ----
# survival, mice, broom, riskRegression, prodlim, scam, mgcv, rms

# Prepared data ----
# data <- readRDS("data/inputs/prepared_cohort.rds")
# Prepared incident cohort; Benzene is the exposure variable.

# Cox sensitivity analyses ----
# cox <- sens_cox(data, "exclude_first_2_years")
sens_cox <- function(data, scenario = c("complete_cases", "exclude_first_2_years",
    "exclude_first_4_years", "exclude_first_8_years", "residence_over_5_years",
    "one_recorded_address", "age_time_scale")) {
  scenario <- match.arg(scenario)
  original <- sens_original(data)
  d <- if (scenario == "complete_cases") {
    original[stats::complete.cases(original[sens_variables()$all]),,drop=FALSE]
  } else sens_impute(original)
  if (grepl("^exclude_first_",scenario)) {
    years <- as.numeric(sub("exclude_first_([248])_years","\\1",scenario))
    keep <- !(d$event==1 & d$time<years)
    d <- d[!is.na(keep)&keep,,drop=FALSE]
  }
  if (scenario == "residence_over_5_years") {
    if (!"p699_i0" %in% names(d)) stop("Supply p699_i0: baseline years at current address.")
    value <- trimws(as.character(d$p699_i0)); years <- rep(NA_real_,length(value))
    numeric_value <- !is.na(value)&grepl("^[0-9]+(\\.[0-9]+)?$",value)
    years[numeric_value] <- as.numeric(value[numeric_value])
    d <- d[!is.na(years)&years>5,,drop=FALSE]
  }
  if (scenario == "one_recorded_address") {
    if (!"n_address_records" %in% names(d)) {
      columns <- grep("^p32220_a[0-9]+$",names(d),value=TRUE)
      if (!length(columns)) stop("Supply n_address_records or p32220_a* address-date columns.")
      d$n_address_records <- rowSums(!is.na(d[columns])&d[columns]!="")
    }
    d <- d[!is.na(d$n_address_records)&d$n_address_records==1,,drop=FALSE]
  }
  d <- sens_prepare(d,impute=FALSE)
  adjustment <- sens_variables()$all
  response <- "survival::Surv(time, event)"
  if (scenario == "age_time_scale") {
    d <- d[is.finite(d$age)&d$age>0,,drop=FALSE]
    d <- sens_prepare(d,impute=FALSE)
    d$age_entry <- d$age; d$age_exit <- d$age+d$time
    adjustment <- setdiff(adjustment,"age")
    response <- "survival::Surv(age_entry, age_exit, event)"
  }
  fits <- lapply(c("exposure_tertile","exposure_SD","exposure_score"),
    function(term) sens_fit(d,term,adjustment,response))
  names(fits) <- c("tertile","sd","trend")
  list(results=sens_model_results(d,fits,scenario),models=fits,
    PH=lapply(fits[c("tertile","sd")],survival::cox.zph),
    tertile_cuts=attr(d,"tertile_cuts"), excluded=nrow(data)-nrow(d))
}

# Multiple imputation ----
# mice <- sens_mice(data)
sens_mice <- function(data, m = 5L, maxit = 10L, seed = 20260910L) {
  d <- sens_prepare(sens_original(data),impute=FALSE)
  variables <- sens_variables()
  imp_data <- d[c(variables$all,"Benzene","event","time")]
  imp_data$nelson_aalen <- mice::nelsonaalen(d,timevar="time",statusvar="event")
  initial <- mice::mice(imp_data,maxit=0,printFlag=FALSE)
  method <- initial$method
  for (v in variables$all) {
    method[v] <- if (!anyNA(imp_data[[v]])) "" else if (v %in% variables$continuous) "pmm" else
      if (nlevels(droplevels(imp_data[[v]]))==2L) "logreg" else "polyreg"
  }
  method[setdiff(names(method),variables$all)] <- ""
  predictor <- initial$predictorMatrix; predictor[,"time"] <- 0
  imputed <- mice::mice(imp_data,m=m,maxit=maxit,method=method,predictorMatrix=predictor,
    seed=seed,printFlag=FALSE,nnet.MaxNWts=10000)
  fits <- lapply(c("exposure_tertile","exposure_SD","exposure_score"),function(term) {
    lapply(seq_len(m),function(i) {
      completed <- mice::complete(imputed,i)
      completed$exposure_tertile <- d$exposure_tertile
      completed$exposure_SD <- d$exposure_SD; completed$exposure_score <- d$exposure_score
      sens_fit(completed,term)
    })
  })
  names(fits) <- c("tertile","sd","trend")
  pooled <- lapply(fits,function(models) summary(mice::pool(mice::as.mira(models))))
  list(results=sens_model_results(d,pooled,"Multiple imputation",pooled=TRUE),
    pooled=pooled,models=fits,imputations=imputed,method=method,predictor_matrix=predictor,
    diagnostics=list(chain_mean=imputed$chainMean,chain_variance=imputed$chainVar,
      logged_events=imputed$loggedEvents),tertile_cuts=attr(d,"tertile_cuts"),
    sample=sens_counts(d))
}

# Competing risks ----
# finegray <- sens_finegray(data)
sens_finegray <- function(data) {
  d <- sens_prepare(sens_impute(sens_original(data)),impute=FALSE)
  death <- as.Date(d$death_date); end <- as.Date(d$All_CB_end_date)
  stopifnot(!anyNA(end))
  d$status_cr <- ifelse(d$event==1,1L,ifelse(!is.na(death)&death<=end,2L,0L))
  Hist <- prodlim::Hist
  fits <- lapply(c("exposure_tertile","exposure_SD","exposure_score"),function(term) {
    f <- sens_formula(term,response="Hist(time, status_cr)"); environment(f) <- environment()
    model <- riskRegression::FGR(f,data=d,cause=1)
    stopifnot(isTRUE(model$crrFit$converged)); model
  })
  names(fits) <- c("tertile","sd","trend")
  tables <- lapply(fits,function(fit) {
    model <- fit$crrFit
    data.frame(term=names(model$coef),estimate=unname(model$coef),
      std.error=sqrt(diag(model$var)),df=Inf)
  })
  list(results=sens_model_results(d,tables,"Fine-Gray",effect="sHR",pooled=TRUE),
    models=fits,status_counts=table(factor(d$status_cr,levels=0:2)),
    tertile_cuts=attr(d,"tertile_cuts"))
}

# Shape-constrained spline ----
# scam <- sens_scam(data)
sens_scam <- function(data) {
  d <- sens_impute(sens_original(data)); d$x <- d$Benzene
  limits <- stats::quantile(d$x[is.finite(d$x)],c(.025,.975),type=7)
  d <- droplevels(d[is.finite(d$x)&d$x>=limits[1]&d$x<=limits[2]&
    is.finite(d$time)&d$time>0&!is.na(d$event),,drop=FALSE])
  stopifnot(all(d$event%in%0:1),sum(d$event)>0,!anyNA(d[sens_variables()$all]))
  reference <- stats::median(d$x); knots <- unname(stats::quantile(d$x,c(.1,.5,.9),type=7))
  stopifnot(length(unique(knots))==3L)
  cuts <- seq(0,ceiling(max(d$time)),by=1); cuts <- cuts[cuts>0&cuts<max(d$time)]
  Surv <- survival::Surv
  ped <- survival::survSplit(Surv(time,event)~.,
    data=d[c("eid","time","event","x",sens_variables()$all)],
    cut=cuts,start="tstart",episode="interval")
  ped$duration <- ped$time-ped$tstart; ped$interval <- factor(ped$interval)
  stopifnot(all(ped$duration>0),sum(ped$event)==sum(d$event),
            isTRUE(all.equal(sum(ped$duration),sum(d$time))))
  s <- mgcv::s; rcs <- rms::rcs
  f <- stats::as.formula(paste("event ~ interval + s(x, bs='mpi', k=6) +",
    paste(sens_variables()$all,collapse=" + "),"+ offset(log(duration))"))
  scam_model <- scam::scam(f,data=ped,family=stats::poisson(link="log"))
  stopifnot(isTRUE(scam_model$conv),all(is.finite(stats::coef(scam_model))))
  f_rcs <- sens_formula(paste0("rcs(x,c(",paste(sprintf("%.17g",knots),collapse=","),"))"))
  environment(f_rcs) <- environment()
  rcs_model <- rms::cph(f_rcs,data=d,x=TRUE,y=TRUE,surv=TRUE)
  stopifnot(!isTRUE(rcs_model$fail),all(is.finite(stats::coef(rcs_model))))
  grid <- sort(unique(c(seq(min(d$x),max(d$x),length.out=200),reference)))
  newdata <- ped[rep(1,length(grid)),,drop=FALSE]; newdata$x <- grid; newdata$duration <- 1
  for (v in sens_variables()$continuous) newdata[[v]] <- stats::median(d[[v]])
  i <- match(reference,grid)
  x_scam <- stats::predict(scam_model,newdata,type="lpmatrix")
  contrast_scam <- sweep(x_scam,2,x_scam[i,],"-")
  eta <- as.numeric(stats::predict(scam_model,newdata,type="link")); beta_scam <- eta-eta[i]
  stopifnot(ncol(x_scam)==nrow(scam_model$Vp.t))
  se_scam <- sqrt(pmax(0,rowSums((contrast_scam%*%scam_model$Vp.t)*contrast_scam)))
  x_rcs <- stats::predict(rcs_model,newdata,type="x")
  contrast_rcs <- sweep(x_rcs,2,x_rcs[i,],"-")
  beta_rcs <- drop(contrast_rcs%*%stats::coef(rcs_model))
  se_rcs <- sqrt(pmax(0,rowSums((contrast_rcs%*%stats::vcov(rcs_model))*contrast_rcs)))
  stopifnot(all(diff(beta_scam)>=-1e-7),abs(beta_scam[i])<1e-8,abs(beta_rcs[i])<1e-8)
  curve <- rbind(data.frame(model="SCAM",concentration=grid,HR=exp(beta_scam),
    lower=exp(beta_scam-1.96*se_scam),upper=exp(beta_scam+1.96*se_scam)),
    data.frame(model="RCS",concentration=grid,HR=exp(beta_rcs),
    lower=exp(beta_rcs-1.96*se_rcs),upper=exp(beta_rcs+1.96*se_rcs)))
  list(curve=curve,models=list(SCAM=scam_model,RCS=rcs_model),sample=sens_counts(d),
    reference=reference,trim_limits=limits,knots=knots,expanded_intervals=nrow(ped))
}

# Helper functions ----
sens_variables <- function() {
  continuous <- c("age", "BMI", "TDI", "sedentary_time", "passive_smoke_home",
                  "outdoor_summer", "outdoor_winter")
  categorical <- c("sex", "ethnicity_2g", "education_2g", "drink_status_2g",
    "smoke_status_2g", "PA_4g", "diet_three", "baseline_hypertension",
    "baseline_dyslipidemia", "baseline_diabetes", "employment_3g")
  list(continuous = continuous, categorical = categorical,
       all = c(continuous, categorical))
}

sens_impute <- function(data) {
  variables <- sens_variables()
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

sens_tertiles <- function(x, cuts = NULL) {
  if (is.null(cuts)) cuts <- unname(stats::quantile(x, c(1/3, 2/3), type = 7))
  stopifnot(length(cuts) == 2L, all(is.finite(cuts)), cuts[1] < cuts[2])
  list(group = cut(x, c(-Inf, cuts, Inf), labels = c("T1", "T2", "T3"),
                   right = TRUE), cuts = cuts)
}

sens_prepare <- function(data, outcome = "All_CB", impute = TRUE) {
  if (impute) data <- sens_impute(data)
  data$exposure <- data$Benzene
  data$time <- data[[paste0(outcome, "_followup_days")]] / 365.25
  data$event <- data[[outcome]]
  stopifnot(is.numeric(data$event), all(is.na(data$event) | data$event %in% 0:1))
  data <- droplevels(data[is.finite(data$exposure) & is.finite(data$time) &
    data$time > 0 & !is.na(data$event), , drop = FALSE])
  stopifnot(nrow(data) > 1L, sum(data$event) > 0)
  tertiles <- sens_tertiles(data$exposure)
  exposure_sd <- stats::sd(data$exposure)
  stopifnot(is.finite(exposure_sd), exposure_sd > 0)
  data$exposure_tertile <- tertiles$group
  data$exposure_score <- as.integer(data$exposure_tertile)
  data$exposure_SD <- data$exposure / exposure_sd
  attr(data, "tertile_cuts") <- tertiles$cuts
  attr(data, "exposure_sd") <- exposure_sd
  data
}

sens_formula <- function(term, covariates = sens_variables()$all,
                           response = "survival::Surv(time, event)") {
  stats::as.formula(paste(response, "~", paste(c(term, covariates), collapse = " + ")))
}

sens_effect <- function(beta, se, df = Inf) {
  critical <- if (is.finite(df)) stats::qt(.975, df) else stats::qnorm(.975)
  log_p <- if (is.finite(df)) log(2) + stats::pt(abs(beta/se), df,
    lower.tail = FALSE, log.p = TRUE) else stats::pchisq((beta/se)^2, 1,
    lower.tail = FALSE, log.p = TRUE)
  data.frame(estimate = unname(beta), SE = unname(se), df = df, HR = exp(unname(beta)),
    lower = exp(unname(beta) - critical * unname(se)),
    upper = exp(unname(beta) + critical * unname(se)), P = exp(log_p), log_P = log_p)
}

sens_extract <- function(model, terms = names(stats::coef(model))) {
  beta <- stats::coef(model)[terms]
  se <- sqrt(diag(stats::vcov(model)))[terms]
  cbind(term = terms, sens_effect(beta, se), row.names = NULL)
}

sens_fit <- function(data, term, covariates = sens_variables()$all,
                      response = "survival::Surv(time, event)") {
  model <- survival::coxph(sens_formula(term, covariates, response), data = data,
    ties = "efron", x = TRUE, y = TRUE, na.action = stats::na.fail)
  stopifnot(model$n == nrow(data), all(is.finite(stats::coef(model))))
  model
}

sens_counts <- function(data) {
  data.frame(participants = nrow(data), cases = sum(data$event),
    person_years = sum(data$time), median_followup_years = stats::median(data$time),
    incidence_10000_py = sum(data$event)/sum(data$time)*10000)
}

sens_original <- function(data) {
  for (v in sens_variables()$categorical) {
    z <- as.character(data[[v]])
    z[is.na(z) | z %in% c("", "99", "Missing")] <- NA_character_
    data[[v]] <- factor(z)
  }
  data$time <- data$All_CB_followup_days/365.25; data$event <- data$All_CB
  data
}

sens_model_results <- function(data, fits, label, effect = "HR", pooled = FALSE) {
  rows <- list()
  trend <- if (pooled) fits$trend[fits$trend$term == "exposure_score", , drop=FALSE] else
    sens_extract(fits$trend, "exposure_score")
  trend_p <- if (pooled) sens_effect(trend$estimate, trend$std.error, trend$df)$P else trend$P
  for (group in c("T1", "T2", "T3", "Per SD increase")) {
    key <- if (group == "Per SD increase") "sd" else "tertile"
    term <- if (key == "sd") "exposure_SD" else paste0("exposure_tertile",group)
    if (group == "T1") {
      result <- data.frame(estimate=0,SE=NA_real_,df=Inf,HR=1,lower=1,upper=1,P=NA_real_,log_P=NA_real_)
    } else if (pooled) {
      z <- fits[[key]][fits[[key]]$term == term,,drop=FALSE]
      stopifnot(nrow(z)==1L)
      result <- sens_effect(z$estimate,z$std.error,z$df)
    } else result <- sens_extract(fits[[key]],term)[,setdiff(names(sens_extract(fits[[key]],term)),"term")]
    selected <- if (group == "Per SD increase") data else data[data$exposure_tertile==group,,drop=FALSE]
    rows[[group]] <- cbind(analysis=label,contrast=group,effect=effect,result,
      sens_counts(selected),SD=stats::sd(data$exposure),
      P_trend=if(group=="Per SD increase") NA_real_ else trend_p)
  }
  do.call(rbind,rows)
}
