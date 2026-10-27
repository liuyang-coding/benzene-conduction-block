# Author: Yang Liu
# Date: 2026-10-27
# R version: 4.5.2

# Packages ----
# survival, rms, interactionR

# Prepared data ----
# data <- readRDS("data/inputs/prepared_cohort.rds")
# Prepared incident cohort; Benzene is the exposure variable.

# PRS associations ----
# prs <- genetic_association(data)
genetic_association <- function(data, spline = TRUE) {
  d <- genetic_prepare(data)
  prs_mean <- mean(d$PRS); prs_sd <- stats::sd(d$PRS)
  stopifnot(prs_sd>0,all(table(d$PRS_tertile)>0))
  d$PRS_z <- (d$PRS-prs_mean)/prs_sd
  models <- list(tertile=genetic_model(d,"PRS_tertile"),per_sd=genetic_model(d,"PRS_z"))
  results <- rbind(cbind(contrast=c("T2 vs T1","T3 vs T1"),
    genetic_extract(models$tertile,c("PRS_tertileT2","PRS_tertileT3"))),
    cbind(contrast="Per SD increase",genetic_extract(models$per_sd,"PRS_z")))
  out <- list(results=results,models=models,PH=lapply(models,survival::cox.zph),
    sample=genetic_counts(d),prs_cuts=attr(d,"prs_cuts"),prs_mean=prs_mean,prs_sd=prs_sd,
    omitted_constants=genetic_covariates(d)$omitted_constants)
  if(spline) {
    limits <- unname(stats::quantile(d$PRS,c(.025,.975),type=7))
    z <- droplevels(d[d$PRS>=limits[1]&d$PRS<=limits[2],,drop=FALSE]); z$prs_raw <- z$PRS
    knots <- unname(stats::quantile(z$prs_raw,c(.10,.50,.90),type=7)); reference <- stats::median(z$prs_raw)
    stopifnot(length(unique(knots))==3L)
    rcs <- rms::rcs
    formula <- genetic_formula(paste0("rcs(prs_raw,c(",paste(sprintf("%.17g",knots),collapse=","),"))"),
      genetic_covariates(z)$active); environment(formula) <- environment()
    fit <- rms::cph(formula,data=z,x=TRUE,y=TRUE,surv=TRUE,method="efron",singular.ok=FALSE)
    stopifnot(!isTRUE(fit$fail),all(is.finite(stats::coef(fit))))
    grid <- sort(unique(c(seq(min(z$prs_raw),max(z$prs_raw),length.out=200),reference)))
    nd <- z[rep(1,length(grid)),,drop=FALSE]; nd$prs_raw <- grid
    design <- stats::predict(fit,nd,type="x"); contrast <- sweep(design,2,design[match(reference,grid),],"-")
    beta <- drop(contrast%*%stats::coef(fit))
    se <- sqrt(pmax(0,rowSums((contrast%*%stats::vcov(fit))*contrast)))
    out$spline <- list(model=fit,anova=stats::anova(fit),PH=survival::cox.zph(fit),
      curve=data.frame(PRS_raw=grid,PRS_z=(grid-prs_mean)/prs_sd,HR=exp(beta),
        lower=exp(beta-1.96*se),upper=exp(beta+1.96*se)),
      reference=reference,knots=knots,trim_limits=limits,sample=genetic_counts(z))
  }
  out
}

# Joint associations ----
# joint <- genetic_joint(data)
genetic_joint <- function(data) {
  d <- genetic_prepare(data,benzene=TRUE)
  joint_levels <- c("T1_T1","T2_T1","T3_T1","T1_T2","T2_T2","T3_T2","T1_T3","T2_T3","T3_T3")
  d$joint <- factor(paste(d$exposure_tertile,d$PRS_tertile,sep="_"),levels=joint_levels)
  if(any(table(d$joint)==0))stop("Empty joint exposure/PRS group.")
  fit <- genetic_model(d,"joint")
  rows <- lapply(joint_levels,function(group) {
    z <- d[d$joint==group,,drop=FALSE]
    result <- if(group=="T1_T1") data.frame(term="Reference",estimate=0,SE=NA_real_,df=Inf,
      HR=1,lower=1,upper=1,P=NA_real_,log_P=NA_real_) else genetic_extract(fit,paste0("joint",group))
    cbind(joint_group=group,result,genetic_counts(z))
  })
  list(results=do.call(rbind,rows),model=fit,PH=survival::cox.zph(fit),interaction=genetic_lrt(d),
    prs_cuts=attr(d,"prs_cuts"),exposure_cuts=attr(d,"exposure_cuts"),
    omitted_constants=attr(fit,"omitted_constants"))
}

# Multiplicative and additive interactions ----
# interaction <- genetic_interaction(data)
genetic_interaction <- function(data) {
  d <- genetic_prepare(data,benzene=TRUE)
  multiplicative <- genetic_lrt(d); additive <- list(); models <- list()
  for(prs_level in c("T2","T3"))for(exposure_level in c("T2","T3")) {
    z <- droplevels(d[d$PRS_tertile%in%c("T1",prs_level)&
      d$exposure_tertile%in%c("T1",exposure_level),,drop=FALSE])
    z$exposure_binary <- as.integer(z$exposure_tertile==exposure_level)
    z$PRS_binary <- as.integer(z$PRS_tertile==prs_level)
    if(any(table(factor(z$exposure_binary,levels=0:1),factor(z$PRS_binary,levels=0:1))==0))
      stop("Empty binary interaction cell.")
    fit <- genetic_model(z,"exposure_binary * PRS_binary")
    result <- interactionR::interactionR(fit,exposure_names=c("exposure_binary","PRS_binary"),
      ci.type="mover",ci.level=.95,em=FALSE,recode=FALSE)$dframe
    result <- result[result$Measures%in%c("RERI","AP"),,drop=FALSE]
    key <- paste(prs_level,exposure_level,sep="/")
    additive[[key]] <- cbind(PRS_group=prs_level,exposure_group=exposure_level,
      participants=nrow(z),cases=sum(z$event),
      HR10=unname(exp(stats::coef(fit)["exposure_binary"])),
      HR01=unname(exp(stats::coef(fit)["PRS_binary"])),
      HR11=exp(sum(stats::coef(fit)[c("exposure_binary","PRS_binary","exposure_binary:PRS_binary")])),result)
    models[[key]] <- fit
  }
  list(multiplicative=multiplicative,additive=do.call(rbind,additive),models=models,
    prs_cuts=attr(d,"prs_cuts"),exposure_cuts=attr(d,"exposure_cuts"))
}

# Helper functions ----
genetic_variables <- function() {
  continuous <- c("age", "BMI", "TDI", "sedentary_time", "passive_smoke_home",
                  "outdoor_summer", "outdoor_winter")
  categorical <- c("sex", "ethnicity_2g", "education_2g", "drink_status_2g",
    "smoke_status_2g", "PA_4g", "diet_three", "baseline_hypertension",
    "baseline_dyslipidemia", "baseline_diabetes", "employment_3g")
  list(continuous = continuous, categorical = categorical,
       all = c(continuous, categorical))
}

genetic_impute <- function(data) {
  variables <- genetic_variables()
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

genetic_tertiles <- function(x, cuts = NULL) {
  if (is.null(cuts)) cuts <- unname(stats::quantile(x, c(1/3, 2/3), type = 7))
  stopifnot(length(cuts) == 2L, all(is.finite(cuts)), cuts[1] < cuts[2])
  list(group = cut(x, c(-Inf, cuts, Inf), labels = c("T1", "T2", "T3"),
                   right = TRUE), cuts = cuts)
}

genetic_formula <- function(term, covariates = genetic_variables()$all,
                           response = "survival::Surv(time, event)") {
  stats::as.formula(paste(response, "~", paste(c(term, covariates), collapse = " + ")))
}

genetic_effect <- function(beta, se, df = Inf) {
  critical <- if (is.finite(df)) stats::qt(.975, df) else stats::qnorm(.975)
  log_p <- if (is.finite(df)) log(2) + stats::pt(abs(beta/se), df,
    lower.tail = FALSE, log.p = TRUE) else stats::pchisq((beta/se)^2, 1,
    lower.tail = FALSE, log.p = TRUE)
  data.frame(estimate = unname(beta), SE = unname(se), df = df, HR = exp(unname(beta)),
    lower = exp(unname(beta) - critical * unname(se)),
    upper = exp(unname(beta) + critical * unname(se)), P = exp(log_p), log_P = log_p)
}

genetic_extract <- function(model, terms = names(stats::coef(model))) {
  beta <- stats::coef(model)[terms]
  se <- sqrt(diag(stats::vcov(model)))[terms]
  cbind(term = terms, genetic_effect(beta, se), row.names = NULL)
}

genetic_counts <- function(data) {
  data.frame(participants = nrow(data), cases = sum(data$event),
    person_years = sum(data$time), median_followup_years = stats::median(data$time),
    incidence_10000_py = sum(data$event)/sum(data$time)*10000)
}

genetic_prepare <- function(data, benzene = FALSE) {
  pcs <- paste0("p22009_a",1:10)
  stopifnot(is.numeric(data$PRS),all(is.finite(data$PRS)))
  prs <- genetic_tertiles(data$PRS)
  d <- genetic_impute(data); d$PRS_tertile <- prs$group
  d$genetic_batch <- factor(d$genetic_batch)
  if (anyNA(d$genetic_batch)||any(!is.finite(as.matrix(d[pcs]))))
    stop("Genetic batch and PCs must be complete.")
  d$time <- d$All_CB_followup_days/365.25; d$event <- d$All_CB
  keep <- is.finite(d$time)&d$time>0&!is.na(d$event)
  if(benzene) keep <- keep&is.finite(d$Benzene)
  d <- droplevels(d[keep,,drop=FALSE]); stopifnot(all(d$event%in%0:1),sum(d$event)>0)
  if(benzene) {
    exposure <- genetic_tertiles(d$Benzene)
    d$exposure_tertile <- exposure$group
    attr(d,"exposure_cuts") <- exposure$cuts
  }
  attr(d,"prs_cuts") <- prs$cuts
  d
}

genetic_covariates <- function(data) {
  all <- c(genetic_variables()$all,"genetic_batch",paste0("p22009_a",1:10))
  active <- all[vapply(data[all],function(x)length(unique(x))>1L,logical(1))]
  list(active=active,omitted_constants=setdiff(all,active))
}

genetic_model <- function(data,term) {
  covariates <- genetic_covariates(data)
  fit <- survival::coxph(genetic_formula(term,covariates$active),data=data,
    ties="efron",na.action=stats::na.fail,singular.ok=FALSE,x=TRUE,y=TRUE)
  stopifnot(fit$n==nrow(data),all(is.finite(stats::coef(fit))),all(is.finite(stats::vcov(fit))))
  attr(fit,"omitted_constants") <- covariates$omitted_constants
  fit
}

genetic_lrt <- function(data) {
  additive <- genetic_model(data,"exposure_tertile + PRS_tertile")
  interaction <- genetic_model(data,"exposure_tertile * PRS_tertile")
  df <- length(stats::coef(interaction))-length(stats::coef(additive))
  statistic <- 2*(as.numeric(stats::logLik(interaction))-as.numeric(stats::logLik(additive)))
  stopifnot(df==4L,statistic>-1e-6)
  list(result=data.frame(chi_square=max(0,statistic),df=df,
    P=stats::pchisq(max(0,statistic),df,lower.tail=FALSE)),
    models=list(additive=additive,interaction=interaction))
}
