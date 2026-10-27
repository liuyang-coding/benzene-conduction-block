# Author: Yang Liu
# Date: 2026-10-27
# R version: 4.5.2

# Packages ----
# survival, rstpm2, WeightIt

# Prepared data ----
# data <- readRDS("data/inputs/prepared_cohort.rds")
# Prepared incident cohort; Benzene is the exposure variable.

# CB-free life expectancy ----
# le <- additional_life_expectancy(data)
additional_life_expectancy <- function(data, grouped = TRUE, reference_age = 45, max_age = 100) {
  stopifnot(max_age>reference_age)
  d <- data; d$time_years <- d$All_CB_followup_days/365.25
  death <- as.Date(d$death_date); end <- as.Date(d$All_CB_end_date)
  stopifnot(!anyNA(end))
  d$CB_free_event <- as.integer(d$All_CB==1|(!is.na(death)&death<=end))
  d <- d[is.finite(d$time_years)&d$time_years>0&!is.na(d$CB_free_event),,drop=FALSE]
  d <- additional_impute(d); cuts <- NULL
  if(grouped) {
    d <- d[is.finite(d$Benzene),,drop=FALSE]
    tertiles <- additional_tertiles(d$Benzene); cuts <- tertiles$cuts
    d$exposure_3g <- factor(c("Low","Medium","High")[as.integer(tertiles$group)],
                            levels=c("High","Medium","Low"))
  }
  f <- additional_formula(if(grouped)"exposure_3g" else character(),
    response="survival::Surv(time_years, CB_free_event)")
  gsm <- rstpm2::gsm
  nsx <- rstpm2::nsx
  environment(f) <- environment()
  fit <- rstpm2::stpm2(f,data=d,df=4)
  beta <- stats::coef(fit); V <- stats::vcov(fit)
  stopifnot(all(is.finite(beta)),all(is.finite(V)))
  if(!is.null(fit@details$convergence)&&fit@details$convergence!=0)stop("stpm2 did not converge.")
  ages <- seq(reference_age,max_age,by=.1); times <- ages-reference_age
  profile <- d[1,additional_variables()$all,drop=FALSE]
  for(v in additional_variables()$continuous)profile[[v]] <- stats::median(d[[v]])
  for(v in additional_variables()$categorical)
    profile[[v]] <- factor(names(which.max(table(d[[v]]))),levels=levels(d[[v]]))
  profile$age <- reference_age
  groups <- if(grouped)levels(d$exposure_3g) else "Overall"
  summaries <- list(); curves <- list(); gradients <- list()
  for(group in groups) {
    nd <- profile[rep(1,length(times)-1L),,drop=FALSE]
    if(grouped)nd$exposure_3g <- factor(group,levels=groups)
    nd$time_years <- times[-1]
    X <- rstpm2::predict(fit,nd,type="lpmatrix"); eta <- drop(X%*%beta)
    H <- exp(eta); S <- exp(-H); survival <- c(1,S)
    if(any(diff(survival)>1e-7))stop("Extrapolated survival increases; inspect spline tail.")
    stopifnot(max(abs(S-as.numeric(rstpm2::predict(fit,nd,type="surv"))))<1e-7)
    se_eta <- sqrt(pmax(0,rowSums((X%*%V)*X)))
    multiplier <- -H*S; multiplier[is.infinite(H)&S==0] <- 0
    gradient <- apply(rbind(rep(0,length(beta)),X*multiplier),2,
      function(g)additional_integral(times,g))
    le <- additional_integral(times,survival)
    se <- sqrt(max(0,drop(t(gradient)%*%V%*%gradient)))
    observed <- max(d$time_years)
    fraction <- if(observed>=max(times))0 else {
      beyond_times <- c(observed,times[times>observed])
      beyond_survival <- c(stats::approx(times,survival,xout=observed)$y,survival[times>observed])
      additional_integral(beyond_times,beyond_survival)/le
    }
    summaries[[group]] <- data.frame(group=group,reference_age=reference_age,max_age=max_age,
      CB_free_years=le,SE=se,lower=le-1.96*se,upper=le+1.96*se,
      maximum_observed_followup=observed,fraction_area_beyond_observed=fraction)
    curves[[group]] <- data.frame(group=group,age=ages,survival=survival,
      lower=c(1,exp(-exp(eta+1.96*se_eta))),upper=c(1,exp(-exp(eta-1.96*se_eta))))
    gradients[[group]] <- gradient
  }
  results <- do.call(rbind,summaries)
  if(grouped)for(group in groups) {
    gradient <- gradients[[group]]-gradients$High
    se <- sqrt(max(0,drop(t(gradient)%*%V%*%gradient)))
    difference <- summaries[[group]]$CB_free_years-summaries$High$CB_free_years
    i <- match(group,results$group)
    results$years_gained[i] <- difference; results$gain_SE[i] <- se
    results$gain_lower[i] <- difference-1.96*se; results$gain_upper[i] <- difference+1.96*se
  }
  list(results=results,survival_curves=do.call(rbind,curves),model=fit,profile=profile,
    tertile_cuts=cuts,participants=nrow(d),first_CB_or_death=sum(d$CB_free_event))
}

# Population attributable fraction ----
# paf <- additional_paf(data)
additional_paf <- function(data) {
  d <- additional_analysis(data)
  fit <- additional_fit(d,"exposure_tertile")
  terms <- c("exposure_tertileT2","exposure_tertileT3")
  b <- stats::coef(fit)[terms]; V <- stats::vcov(fit)[terms,terms,drop=FALSE]
  hr <- exp(b); prevalence <- as.numeric(prop.table(table(d$exposure_tertile)))[2:3]
  prevalence_cov <- (diag(prevalence)-tcrossprod(prevalence))/nrow(d)
  rows <- lapply(list(T2=1L,T3=2L,T2_and_T3=1:2),function(index) {
    p <- prevalence[index]; h <- hr[index]; denominator <- 1+sum(p*(h-1))
    value <- 1-1/denominator
    gradient_b <- p*h/denominator^2; gradient_p <- (h-1)/denominator^2
    variance <- drop(t(gradient_b)%*%V[index,index,drop=FALSE]%*%gradient_b)+
      drop(t(gradient_p)%*%prevalence_cov[index,index,drop=FALSE]%*%gradient_p)
    se <- sqrt(max(0,variance))
    data.frame(prevalence=sum(p),PAF_percent=100*value,SE_percent=100*se,
      lower_percent=100*(value-1.96*se),upper_percent=100*(value+1.96*se))
  })
  results <- do.call(rbind,rows); results$contrast <- rownames(results)
  list(results=results,HR=additional_extract(fit,terms),model=fit,PH=survival::cox.zph(fit),
    sample=additional_counts(d),tertile_cuts=attr(d,"tertile_cuts"))
}

# Potential impact fraction ----
# pif <- additional_pif(data)
additional_pif <- function(data, n_grid = 180L) {
  d <- additional_analysis(data)
  fit <- additional_fit(d,"exposure")
  beta <- unname(stats::coef(fit)["exposure"]); se <- sqrt(stats::vcov(fit)["exposure","exposure"])
  x <- d$exposure; probs <- c(.75,.50,.25,.10)
  targets <- unname(stats::quantile(x,probs,type=7)); stopifnot(min(targets)>0,min(targets)<max(x))
  value <- function(target,b) {
    counterfactual <- pmin(x,target); center <- max(c(b*x,b*counterfactual))
    100*(1-mean(exp(b*counterfactual-center))/mean(exp(b*x-center)))
  }
  calculate <- function(target) {
    edges <- c(value(target,beta-1.96*se),value(target,beta+1.96*se))
    data.frame(target_concentration=target,control_intensity_percent=100*(mean(x)-mean(pmin(x,target)))/mean(x),
      population_above_target_percent=100*mean(x>target),PIF_percent=value(target,beta),
      lower_percent=min(edges),upper_percent=max(edges))
  }
  scenarios <- do.call(rbind,lapply(targets,calculate)); scenarios$target_percentile <- 100*probs
  list(scenarios=scenarios,curve=do.call(rbind,lapply(seq(max(x),min(targets),length.out=n_grid),calculate)),
    model=fit,beta=beta,beta_SE=se,SD=stats::sd(x),HR_per_SD=exp(beta*stats::sd(x)),
    sample=additional_counts(d))
}

# Generalized propensity scores ----
# gps <- additional_gps(data)
additional_gps <- function(data, seed = 20260910L) {
  d <- additional_analysis(data)
  covariates <- additional_variables()$all
  f_gps <- stats::reformulate(covariates,response="exposure")
  f_ipw <- stats::reformulate(covariates,response="exposure_tertile")
  set.seed(seed); ipw <- WeightIt::weightit(f_ipw,data=d,method="glm",estimand="ATE",stabilize=TRUE)
  set.seed(seed); gps <- WeightIt::weightit(f_gps,data=d,method="glm")
  for(w in list(ipw$weights,gps$weights))stopifnot(length(w)==nrow(d),all(is.finite(w)),all(w>0))
  density_model <- stats::lm(f_gps,data=d)
  residual_sd <- stats::sigma(density_model); stopifnot(is.finite(residual_sd),residual_sd>0)
  d$GPS <- stats::dnorm(d$exposure,mean=stats::predict(density_model),sd=residual_sd)
  propensity <- ipw$ps
  if(is.null(propensity)) propensity <- as.numeric(prop.table(table(d$exposure_tertile))[as.character(d$exposure_tertile)])/ipw$weights
  if(is.matrix(propensity)||is.data.frame(propensity)) {
    propensity <- as.matrix(propensity)
    index <- match(as.character(d$exposure_tertile),colnames(propensity)); stopifnot(!anyNA(index))
    propensity <- propensity[cbind(seq_len(nrow(d)),index)]
  }
  stopifnot(length(propensity)==nrow(d),all(is.finite(propensity)),all(propensity>0&propensity<=1))
  matrix <- stats::model.matrix(stats::reformulate(covariates,intercept=FALSE),d)
  matrix <- matrix[,apply(matrix,2,function(x)length(unique(x))>1),drop=FALSE]
  lower <- c(0,.005,.01,.025); upper <- 1-lower
  labels <- c("No truncation","0.5-99.5 percentiles","1-99 percentiles","2.5-97.5 percentiles")
  results <- list(); models <- list(); balance <- list(); weights_summary <- list()
  for(i in seq_along(lower)) {
    clipped <- lapply(list(IPW=ipw$weights,GPS=gps$weights),function(w) {
      limits <- stats::quantile(w,c(lower[i],upper[i]),names=FALSE,type=7)
      pmin(pmax(w,limits[1]),limits[2])
    })
    for(framework in names(clipped)) {
      w <- clipped[[framework]]
      b <- additional_balance(matrix,d$exposure,d$exposure_tertile,w,framework)
      b$unweighted <- additional_balance(matrix,d$exposure,d$exposure_tertile,rep(1,nrow(d)),framework)$balance
      b$scenario <- labels[i]; b$below_0.10 <- b$balance<.10
      balance[[paste(i,framework)]] <- b
      weights_summary[[paste(i,framework)]] <- data.frame(scenario=labels[i],framework=framework,
        min=min(w),median=stats::median(w),max=max(w),ESS=sum(w)^2/sum(w^2),
        ESS_percent=100*sum(w)^2/sum(w^2)/nrow(d))
    }
    d$analysis_weight <- clipped$GPS
    fit <- survival::coxph(additional_formula("exposure_SD"),data=d,weights=analysis_weight,
      robust=TRUE,cluster=eid,x=TRUE,y=TRUE,model=TRUE,na.action=stats::na.fail,ties="efron")
    stopifnot(fit$n==nrow(d),all(is.finite(stats::coef(fit))))
    results[[labels[i]]] <- cbind(method="GPS-weighted Cox, 18 covariates",scenario=labels[i],
      additional_extract(fit,"exposure_SD")); models[[labels[i]]] <- fit
  }
  fit <- survival::coxph(additional_formula("exposure_SD + GPS",character()),data=d,
    robust=TRUE,cluster=eid,x=TRUE,y=TRUE,model=TRUE,na.action=stats::na.fail,ties="efron")
  stopifnot(fit$n==nrow(d),all(is.finite(stats::coef(fit))))
  results$linear_GPS <- cbind(method="Unweighted Cox with linear GPS",scenario="Not applicable",
    additional_extract(fit,"exposure_SD")); models$linear_GPS <- fit
  level <- ifelse(d$exposure<=stats::median(d$exposure),"Low","High")
  ranges <- t(vapply(split(d$GPS,level),range,numeric(2)))
  common <- c(max(ranges[,1]),min(ranges[,2]))
  support <- data.frame(eid=d$eid,GPS=d$GPS,exposure_level=level,
    observed_treatment_probability=as.numeric(propensity),
    outside_common_support=d$GPS<common[1]|d$GPS>common[2])
  list(results=do.call(rbind,results),models=models,PH=lapply(models,survival::cox.zph),
    weights_summary=do.call(rbind,weights_summary),balance=do.call(rbind,balance),
    propensity=list(IPW=ipw,GPS=gps),support=support,common_support=common,
    sample=additional_counts(d),SD=attr(d,"exposure_sd"),tertile_cuts=attr(d,"tertile_cuts"))
}

# Helper functions ----
additional_variables <- function() {
  continuous <- c("age", "BMI", "TDI", "sedentary_time", "passive_smoke_home",
                  "outdoor_summer", "outdoor_winter")
  categorical <- c("sex", "ethnicity_2g", "education_2g", "drink_status_2g",
    "smoke_status_2g", "PA_4g", "diet_three", "baseline_hypertension",
    "baseline_dyslipidemia", "baseline_diabetes", "employment_3g")
  list(continuous = continuous, categorical = categorical,
       all = c(continuous, categorical))
}

additional_impute <- function(data) {
  variables <- additional_variables()
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

additional_tertiles <- function(x, cuts = NULL) {
  if (is.null(cuts)) cuts <- unname(stats::quantile(x, c(1/3, 2/3), type = 7))
  stopifnot(length(cuts) == 2L, all(is.finite(cuts)), cuts[1] < cuts[2])
  list(group = cut(x, c(-Inf, cuts, Inf), labels = c("T1", "T2", "T3"),
                   right = TRUE), cuts = cuts)
}

additional_prepare <- function(data, outcome = "All_CB", impute = TRUE) {
  if (impute) data <- additional_impute(data)
  data$exposure <- data$Benzene
  data$time <- data[[paste0(outcome, "_followup_days")]] / 365.25
  data$event <- data[[outcome]]
  stopifnot(is.numeric(data$event), all(is.na(data$event) | data$event %in% 0:1))
  data <- droplevels(data[is.finite(data$exposure) & is.finite(data$time) &
    data$time > 0 & !is.na(data$event), , drop = FALSE])
  stopifnot(nrow(data) > 1L, sum(data$event) > 0)
  tertiles <- additional_tertiles(data$exposure)
  exposure_sd <- stats::sd(data$exposure)
  stopifnot(is.finite(exposure_sd), exposure_sd > 0)
  data$exposure_tertile <- tertiles$group
  data$exposure_score <- as.integer(data$exposure_tertile)
  data$exposure_SD <- data$exposure / exposure_sd
  attr(data, "tertile_cuts") <- tertiles$cuts
  attr(data, "exposure_sd") <- exposure_sd
  data
}

additional_formula <- function(term, covariates = additional_variables()$all,
                           response = "survival::Surv(time, event)") {
  stats::as.formula(paste(response, "~", paste(c(term, covariates), collapse = " + ")))
}

additional_effect <- function(beta, se, df = Inf) {
  critical <- if (is.finite(df)) stats::qt(.975, df) else stats::qnorm(.975)
  log_p <- if (is.finite(df)) log(2) + stats::pt(abs(beta/se), df,
    lower.tail = FALSE, log.p = TRUE) else stats::pchisq((beta/se)^2, 1,
    lower.tail = FALSE, log.p = TRUE)
  data.frame(estimate = unname(beta), SE = unname(se), df = df, HR = exp(unname(beta)),
    lower = exp(unname(beta) - critical * unname(se)),
    upper = exp(unname(beta) + critical * unname(se)), P = exp(log_p), log_P = log_p)
}

additional_extract <- function(model, terms = names(stats::coef(model))) {
  beta <- stats::coef(model)[terms]
  se <- sqrt(diag(stats::vcov(model)))[terms]
  cbind(term = terms, additional_effect(beta, se), row.names = NULL)
}

additional_fit <- function(data, term, covariates = additional_variables()$all,
                      response = "survival::Surv(time, event)") {
  model <- survival::coxph(additional_formula(term, covariates, response), data = data,
    ties = "efron", x = TRUE, y = TRUE, na.action = stats::na.fail)
  stopifnot(model$n == nrow(data), all(is.finite(stats::coef(model))))
  model
}

additional_counts <- function(data) {
  data.frame(participants = nrow(data), cases = sum(data$event),
    person_years = sum(data$time), median_followup_years = stats::median(data$time),
    incidence_10000_py = sum(data$event)/sum(data$time)*10000)
}

additional_analysis <- function(data) {
  time <- data$All_CB_followup_days/365.25
  data <- data[is.finite(time)&time>0&!is.na(data$All_CB),,drop=FALSE]
  additional_prepare(additional_impute(data),impute=FALSE)
}

additional_integral <- function(x,y) sum(diff(x)*(head(y,-1)+tail(y,-1)))/2

additional_weighted_correlation <- function(x,y,w) {
  mx <- sum(w*x)/sum(w); my <- sum(w*y)/sum(w)
  vx <- sum(w*(x-mx)^2)/sum(w); vy <- sum(w*(y-my)^2)/sum(w)
  if(vx<=0||vy<=0)return(NA_real_)
  (sum(w*(x-mx)*(y-my))/sum(w))/sqrt(vx*vy)
}

additional_balance <- function(matrix,exposure,group,w,framework) {
  values <- vapply(seq_len(ncol(matrix)),function(j) {
    x <- matrix[,j]
    if(framework=="GPS")return(abs(additional_weighted_correlation(exposure,x,w)))
    pairs <- utils::combn(levels(group),2)
    max(apply(pairs,2,function(pair) {
      a <- group==pair[1]; b <- group==pair[2]
      difference <- sum(w[a]*x[a])/sum(w[a])-sum(w[b]*x[b])/sum(w[b])
      denominator <- sqrt((mean((x[a]-mean(x[a]))^2)+mean((x[b]-mean(x[b]))^2))/2)
      if(denominator<=0) { if(abs(difference)<1e-12)0 else Inf } else abs(difference)/denominator
    }))
  },numeric(1))
  data.frame(covariate=colnames(matrix),balance=values,framework=framework)
}
