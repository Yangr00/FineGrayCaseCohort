#' Fit and diagnose a case-cohort subdistribution hazards model
#'
#' Fits the proportional subdistribution hazards model once, then computes three
#' proportionality diagnostics. All covariates enter through the same formula.
#' The input must retain every member of the full cohort. Measurements are
#' required for subcohort members and target-event cases; covariates may be
#' missing for other subjects. Rows are never silently omitted.
#'
#' @param formula One-sided covariate formula, e.g. `~ covariate1 + covariate2`.
#'   Factors use standard model-matrix contrasts. The intercept is removed.
#'   A dot selects all columns except the specified time, status and subcohort
#'   columns. Use `~ . - id` to exclude an additional non-predictor column.
#' @param data Full-cohort data frame, one row per subject.
#' @param time,status,subcohort Character strings naming columns in `data`.
#'   Time is finite and non-negative; status uses 0 for censoring and positive
#'   integers for event causes. Subcohort is 0/1 or logical random-subcohort
#'   membership, distinct from measurement availability.
#' @param cause Positive integer code for the event of interest. Default 1.
#' @param time_transform Function defining q(t) in the regression alternative
#'   beta(t)=beta+diag(q(t))*theta. Default q(t)=t. Returns a finite numeric vector
#'   of length n or an n-by-p matrix, with p model-matrix columns. A constant
#'   transform has no identifiable time trend and yields an unavailable reg test.
#' @details Only random unstratified subcohort sampling with the implemented
#'   canonical weights is supported. Delayed entry, stratified sampling and
#'   time-varying input covariates are not supported.
#'
#'   Cor is Pearson correlation of event time and residual. Rank cor is Pearson
#'   correlation of ranked event time and residual, with average ranks for ties;
#'   it is not Spearman correlation. Their statistics use a t reference with
#'   n_events-2 degrees of freedom. Reg tests theta=0, with beta estimated under
#'   the null, using a chi-squared reference. With multiple coefficients reg
#'   provides a joint GLOBAL row and individual component rows. Factor columns
#'   are reported separately, not as a joint factor-term test.
#'
#'   Beta is fitted by one call to stats::nlminb with analytic derivatives and
#'   independent convergence checks. In standardized coefficient coordinates,
#'   the maximum absolute score per target event must be at most 1e-7 and the
#'   maximum relative Newton correction at most 1e-6. The final likelihood and
#'   derivatives must be finite, with identifiable positive information. A
#'   nonzero optimizer return code can be accepted when these checks pass;
#'   its original code and message remain in `solver_diagnostics`.
#'   Failure gives a warning and NA inference
#'   together with diagnostic status. Invalid input gives an error. If beta
#'   converges but its covariance is unavailable, the coefficient estimate and
#'   descriptive residuals are retained; covariance matrices and all diagnostic
#'   p-values are NA. The convergence flag refers to beta optimization only.
#' @return An S3 object of class `cc_fit` with components:
#'   \item{coefficients}{Named fitted beta vector.}
#'   \item{var, naive.var}{Sandwich and information covariance matrices.}
#'   \item{tests}{Data frame with method, covariate, estimate (correlation),
#'     statistic, df, p.value, status and message.}
#'   \item{regression_test}{Nuisance-adjusted score and covariance, component p-values,
#'     statistic, degrees of freedom and status for the regression test.}
#'   \item{residuals, event_residuals}{Full-cohort and target-event residual
#'     matrices. Non-target rows in the full matrix are NA.}
#'   \item{event_times, event_index}{Times and 1-based input row indices matching
#'     the rows of event_residuals. Tied events remain separate rows.}
#'   \item{n, n_events, n_subcohort}{Full-cohort, target-event and subcohort counts.}
#'   \item{converged, solver_method, solver_status, solver_message}{Beta optimization status.}
#'   \item{solver_diagnostics}{Original optimizer code, message, iterations and
#'     evaluations, plus available score, correction and information checks.
#'     `verified` records whether the independent numerical checks passed.}
#'   \item{variance_status, variance_message}{Covariance availability and diagnostic
#'     message; inference is withheld if beta covariance is unavailable.}
#'   \item{formula, terms, xlevels, contrasts, coefficient_terms}{Model metadata;
#'     coefficient_terms maps expanded columns to formula terms.}
#'   \item{cause, time_transform, call}{Target cause, time transformation and call.}
#' @examples
#' dat <- read.csv(system.file("extdata", "example.csv",
#'                            package = "CompetingRiskRPackage"))
#' fit <- cc_fit(~ covariate1, dat, "time", "status", "subcohort")
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#' head(residuals(fit, type = "event"))
#' @export
cc_fit <- function(formula, data, time, status, subcohort, cause = 1L,
                   time_transform = function(t) t) {
  model <- .prepare_data(formula,data,time,status,subcohort,cause,time_transform)
  out <- .fit_model(model)
  out$cause <- cause
  out$formula <- formula
  out$terms <- model$terms
  out$xlevels <- model$xlevels
  out$contrasts <- model$contrasts
  out$coefficient_terms <- model$coefficient_terms
  out$call <- match.call()
  out
}

# All fitting steps share the prepared full cohort and one censoring distribution.
# details is used only by the test reference adapter, never by cc_fit().
.fit_model <- function(model, details=FALSE) {
  .get_beta_optimizer()
  km <- censoring_create(model$time,1L-model$delta)
  beta <- .fit_beta(model,km)
  labels <- colnames(model$z)
  p <- ncol(model$z)
  out <- list(coefficients=stats::setNames(beta$coefficients,labels),
    n=length(model$time),n_events=length(model$event_index),n_subcohort=sum(model$cc),
    event_index=model$event_index,event_times=model$time[model$event_index],
    converged=beta$status=='converged',solver_method='nlminb',
    solver_status=beta$status,solver_message=beta$message,time_transform=model$time_transform)
  out$solver_diagnostics <- beta$solver_diagnostics
  out$variance_status <- 'not_computed'
  out$variance_message <- 'The null fit did not converge.'
  covariance <- list(var=matrix(NA_real_,p,p),naive.var=matrix(NA_real_,p,p))
  diagnostics <- list(residuals=matrix(NA_real_,out$n,p),cor=rep(NA_real_,p),
    cor_p=rep(NA_real_,p),rank_cor=rep(NA_real_,p),rank_cor_p=rep(NA_real_,p))
  if(out$converged) {
    checked <- tryCatch({
      value <- beta_covariance(model$time,model$event,model$z,model$cc,
                               out$coefficients,km,details)
      .validate_beta_covariance(value,p)
      list(value=value,status='ok',message='')
    },error=function(e) list(value=NULL,status='unavailable',message=conditionMessage(e)))
    out$variance_status <- checked$status
    out$variance_message <- checked$message
    if(checked$status=='ok') covariance <- checked$value else
      warning(paste('Beta covariance unavailable; inference withheld:',checked$message),call.=FALSE)
    diagnostics <- residual_diagnostics(model$time,model$event,model$z,model$cc,
                                         out$coefficients,km)
    bad_residual <- colSums(!is.finite(diagnostics$residuals[model$event_index,,drop=FALSE]))>0
    if(any(bad_residual))
      warning(paste('Non-finite event residuals; correlation inference unavailable for',
                    paste(labels[bad_residual],collapse=', ')),call.=FALSE)
    if(out$n_events<=2L)
      warning('At most two target events: correlation p-values are unavailable; model estimates are still returned.')
  } else warning(out$solver_message,call.=FALSE)
  out$var <- covariance$var
  out$naive.var <- covariance$naive.var
  dimnames(out$var) <- dimnames(out$naive.var) <- list(labels,labels)
  out$residuals <- diagnostics$residuals
  colnames(out$residuals) <- labels
  out$event_residuals <- out$residuals[out$event_index,,drop=FALSE]
  out$regression_test <- .regression_test(model,out,km)
  out$tests <- .diagnostic_table(out,diagnostics)
  class(out) <- 'cc_fit'
  if(details) return(list(fit=out,score_covariance=covariance$score_covariance,censoring=km))
  out
}

# Check complete matrices in dimensionless coordinates, preserving unit invariance.
.validate_beta_covariance <- function(value,p) {
  for(name in c('var','naive.var')) {
    v <- value[[name]]
    if(!is.numeric(v) || !is.matrix(v) || !identical(dim(v),c(p,p)) ||
       any(!is.finite(v)) || any(diag(v)<=0))
      stop(paste(name,'is not a finite covariance matrix with positive diagonal.'))
    scaled <- v/outer(sqrt(diag(v)),sqrt(diag(v)))
    # Retain roundoff-level asymmetry in the canonical near-collinear fixture.
    if(max(abs(scaled-t(scaled)))>1e-7)
      stop(paste(name,'is not symmetric.'))
    ev <- eigen((scaled+t(scaled))/2,symmetric=TRUE,only.values=TRUE)$values
    if(any(!is.finite(ev)) || min(ev)<=max(1e-12,max(ev)*1e-12))
      stop(paste(name,'is singular or not positive definite.'))
  }
  invisible(TRUE)
}
