# Regression has its own joint influence covariance, distinct from vcov(fit).
.regression_test <- function(model,fit,censoring) {
  X <- model$time; event <- model$event; Z <- model$z; cc <- model$cc
  p <- ncol(Z); q <- model$q
  reg <- list(statistic=NA_real_,df=p,p_value=NA_real_,component_p_values=rep(NA_real_,p),
              status="fit_nonconverged",message="The null fit did not converge.")
  if(isTRUE(fit$converged) && fit$variance_status!='ok') {
    reg$status <- 'variance_unavailable'
    reg$message <- fit$variance_message
  } else if (isTRUE(fit$converged)) {
    reg <- tryCatch({
      measured <- cc==1 | event==1
      clean <- Z
      clean[!measured,] <- 0
      scales <- apply(clean[measured,,drop=FALSE],2,stats::sd)
      centers <- colMeans(clean[measured,,drop=FALSE])
      z <- sweep(sweep(clean,2,centers,"-"),2,scales,"/")
      z[!measured,] <- 0
      # Center/scale Q on event times. Adding a constant multiple of U_beta or
      # rescaling theta leaves the nuisance-adjusted score test unchanged.
      for(k in seq_len(p)) {
        qm <- mean(q[event==1,k]);qs <- stats::sd(q[event==1,k])
        if (!is.finite(qs) || qs<=0) stop("The time transform has no variation at target-event times.")
        q[,k] <- (q[,k]-qm)/qs
      }
      components <- regression_score_components(as.numeric(X),as.integer(event),z,
                           as.integer(cc),as.numeric(fit$coefficients)*scales,q,censoring)
      B <- components$Omega_theta_beta %*% solve(components$Omega_beta_beta)
      projection <- cbind(-B,diag(p))
      covariance <- projection %*% components$Sigma %*% t(projection)
      covariance <- (covariance+t(covariance))/2
      ev <- eigen(covariance,symmetric=TRUE,only.values=TRUE)$values
      if(any(!is.finite(ev)) || min(ev)<=max(1e-12,max(ev)*1e-12))
        stop("The nuisance-adjusted regression covariance is singular or not positive definite.")
      efficient_score <- as.numeric(projection %*% components$score)
      statistic <- as.numeric(crossprod(efficient_score,solve(covariance,efficient_score))/length(X))
      list(statistic=statistic,df=p,p_value=stats::pchisq(statistic,p,lower.tail=FALSE),
           component_p_values=stats::pchisq(efficient_score^2/(length(X)*diag(covariance)),1,lower.tail=FALSE),
           status="ok",message="H0: theta=0, with beta estimated under H0; case-cohort and censoring variance corrections.",
           score=efficient_score,covariance=covariance)
    },error=function(e) {
      warning(paste("Regression proportionality test unavailable:",conditionMessage(e)),call.=FALSE)
      list(statistic=NA_real_,df=p,p_value=NA_real_,component_p_values=rep(NA_real_,p),
           status="unavailable",message=conditionMessage(e))
    })
  }
  names(reg$component_p_values) <- colnames(Z)
  reg
}

.diagnostic_table <- function(x,diagnostics) {
  labels <- names(x$coefficients)
  correlation_rows <- function(method, r, p) {
    df <- x$n_events-2L
    valid <- x$converged & x$variance_status=='ok' & is.finite(p) & is.finite(r) & abs(r)<=1 & df>0
    p[!valid] <- NA_real_
    statistic <- rep(NA_real_,length(r))
    statistic[valid] <- r[valid]*sqrt(df/pmax(0,1-r[valid]^2))
    data.frame(method=method,covariate=labels,estimate=as.numeric(r),
      statistic=statistic,df=rep(if(df>0) df else NA_real_,length(r)),
      p.value=as.numeric(p),status=ifelse(valid,"ok","unavailable"),
      message=ifelse(valid,"",if(!x$converged) x$solver_message else
        if(x$variance_status!='ok') x$variance_message else
        "Insufficient target events or no variation for correlation inference."),
      row.names=NULL)
  }
  reg <- x$regression_test
  k <- length(labels)
  regrow <- data.frame(method="reg",covariate=if(k==1L) labels else "GLOBAL",
    estimate=NA_real_,statistic=reg$statistic,df=reg$df,p.value=reg$p_value,
    status=reg$status,message=if(reg$status=="ok") "" else reg$message)
  if(k>1L) {
    stat <- if(reg$status=="ok") reg$score^2/(x$n*diag(reg$covariance)) else rep(NA_real_,k)
    regrow <- rbind(regrow,data.frame(method="reg",covariate=labels,estimate=NA_real_,
      statistic=as.numeric(stat),df=1,p.value=as.numeric(reg$component_p_values),
      status=reg$status,message=if(reg$status=="ok") "" else reg$message))
  }
  tests <- rbind(correlation_rows("cor",diagnostics$cor,diagnostics$cor_p),
    correlation_rows("rank cor",diagnostics$rank_cor,diagnostics$rank_cor_p),regrow)
  rownames(tests) <- NULL
  tests
}

