#' Methods for case-cohort model objects
#' @param object,x A `cc_fit` object (or `summary.cc_fit` for its print method).
#' @param type Residual rows: `"full"` (default) preserves all input rows;
#'   `"event"` selects target-event rows in input order. Both return a matrix.
#' @param ... Reserved for future options.
#' @return `coef` returns a named vector; `vcov` returns the sandwich covariance
#'   matrix; `residuals` returns a matrix. `summary` returns a `summary.cc_fit`
#'   object with call, counts, solver and covariance status, coefficients (estimate and sandwich
#'   and information standard errors), tests and regression-test details.
#'   Print methods invisibly return their input. `plot` prints residual plots
#'   and invisibly returns a named list of ggplot objects, one per coefficient.
#' @name cc_fit_methods
NULL

#' @rdname cc_fit_methods
#' @export
coef.cc_fit <- function(object, ...) object$coefficients
#' @rdname cc_fit_methods
#' @export
vcov.cc_fit <- function(object, ...) object$var
#' @rdname cc_fit_methods
#' @export
residuals.cc_fit <- function(object, type=c("full","event"), ...) {
  type <- match.arg(type)
  if(type=="full") object$residuals else object$event_residuals
}
#' @rdname cc_fit_methods
#' @export
summary.cc_fit <- function(object, ...) {
  out <- list(call=object$call,converged=object$converged,
    solver_status=object$solver_status,solver_method=object$solver_method,
    solver_message=object$solver_message,variance_status=object$variance_status,
    solver_diagnostics=object$solver_diagnostics,
    variance_message=object$variance_message,n=object$n,n_events=object$n_events,
    n_subcohort=object$n_subcohort,
    coefficients=data.frame(estimate=object$coefficients,
      SE_sandwich=sqrt(diag(object$var)),SE_information=sqrt(diag(object$naive.var))),
    tests=object$tests,regression_test=object$regression_test)
  class(out) <- "summary.cc_fit"
  out
}
#' @rdname cc_fit_methods
#' @export
print.summary.cc_fit <- function(x, ...) {
  cat("Case-cohort subdistribution hazards model\n")
  cat("Full cohort:",x$n," Target events:",x$n_events," Subcohort:",x$n_subcohort,"\n")
  if(!isTRUE(x$converged)) cat("Fit status:",x$solver_status,"\n",x$solver_message,"\n")
  if(isTRUE(x$converged) && x$variance_status!='ok')
    cat("Covariance status:",x$variance_status,"\n",x$variance_message,"\nInference is unavailable.\n")
  print(x$coefficients)
  cat("\nProportionality diagnostics (cor/rank cor: t; reg: chi-squared)\n")
  print(x$tests[,c("method","covariate","statistic","df","p.value","status")],row.names=FALSE)
  invisible(x)
}
#' @rdname cc_fit_methods
#' @export
print.cc_fit <- function(x, ...) {
  print(summary(x))
  invisible(x)
}
