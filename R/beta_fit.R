# Internal: analytic callbacks share the C++ evaluation cache at each parameter.
# A single stats::nlminb call is made; no fallback restarts or penalty terms.
.get_beta_optimizer <- function() {
  if (!requireNamespace("stats", quietly = TRUE))
    stop("Required optimizer stats::nlminb is unavailable: R's 'stats' package could not be loaded. Check or repair your R installation.", call. = FALSE)
  optimizer <- tryCatch(getExportedValue("stats", "nlminb"),
                        error = function(e) NULL)
  if (!is.function(optimizer))
    stop("Required optimizer stats::nlminb is unavailable. The nlminb function is supplied by R's bundled 'stats' package. Check or repair your R installation.", call. = FALSE)
  optimizer
}

.fit_beta_nlminb <- function(problem, p) {
  optimizer <- .get_beta_optimizer()
  evaluate <- function(beta_scaled) beta_likelihood_eval(problem, beta_scaled)
  optimizer(rep(0, p),
    objective = function(beta_scaled) evaluate(beta_scaled)$value,
    gradient = function(beta_scaled) evaluate(beta_scaled)$gradient,
    hessian = function(beta_scaled) evaluate(beta_scaled)$hessian,
    control = list(iter.max = 1000, eval.max = 2000,
                   rel.tol = 1e-12, x.tol = 1e-10, sing.tol = 1e-12))
}

# R owns the optimizer lifecycle; native code only evaluates/checks the problem.
.fit_beta <- function(model,censoring) {
  p <- ncol(model$z)
  result <- NULL
  fail <- function(status,message)
    list(coefficients=rep(NA_real_,p),status=status,message=message,
      solver_diagnostics=list(optimizer_run=!is.null(result),optimizer_code=result$convergence,
        optimizer_message=result$message,iterations=result$iterations,
        evaluations=result$evaluations,verified=FALSE))
  tryCatch({
    problem <- beta_likelihood_create(model$time,model$event,model$z,model$cc,censoring)
    on.exit(beta_likelihood_release(problem$pointer),add=TRUE)
    if(problem$separated)
      return(fail('no_finite_estimate','No finite unpenalized estimate: a separating risk-set direction was detected.'))
    result <- .fit_beta_nlminb(problem$pointer,p)
    checked <- beta_likelihood_check(problem$pointer,result)
    checked
  },error=function(e) fail('nonconverged',conditionMessage(e)))
}
