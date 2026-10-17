.prepare_data <- function(formula,data,time,status,subcohort,cause,time_transform) {
  if (!inherits(formula, "formula") || length(formula) != 2L)
    stop("Use a one-sided covariate formula, for example ~ covariate1 + covariate2.")
  if (!is.data.frame(data) || !nrow(data)) stop("data must be a non-empty full-cohort data frame.")
  for (column in list(time, status, subcohort))
    if (!is.character(column) || length(column)!=1L || is.na(column) || !column %in% names(data))
      stop("time, status and subcohort must each name one column in data.")
  if (!is.numeric(cause) || length(cause)!=1L || !is.finite(cause) || cause<1 || cause!=trunc(cause))
    stop("cause must be a positive integer.")
  original <- data[[status]]
  if (!is.numeric(original) || any(!is.finite(original)) || any(original<0 | original!=trunc(original)))
    stop("status must contain non-negative integer codes, with 0 denoting censoring.")
  cc <- data[[subcohort]]
  if (is.logical(cc)) cc <- as.integer(cc)
  coded <- ifelse(original==0, 0L, ifelse(original==cause, 1L, 2L))
  # Expand dots over predictors only; cohort design columns are not predictors.
  formula_data <- if ("." %in% all.names(formula))
    data[,setdiff(names(data),c(time,status,subcohort)),drop=FALSE] else data
  mf <- stats::model.frame(formula, data=formula_data, na.action=stats::na.pass,
                           drop.unused.levels=TRUE)
  design <- stats::model.matrix(attr(mf,"terms"), mf)
  assignment <- attr(design,"assign")
  keep <- colnames(design)!="(Intercept)"
  contrasts <- attr(design,"contrasts")
  design <- design[,keep,drop=FALSE]
  if (nrow(design)!=nrow(data)) stop("Model construction removed cohort members.")
  model <- .prepare_numeric(data[[time]],as.integer(original>0),coded,design,cc,time_transform)
  model$terms <- attr(mf,'terms')
  model$xlevels <- lapply(mf[vapply(mf,is.factor,logical(1))],levels)
  model$contrasts <- contrasts
  term_names <- attr(model$terms,'term.labels')
  model$coefficient_terms <- stats::setNames(term_names[assignment[keep]],colnames(design))
  model
}

# Matrix-level preparation also permits independent teacher-code reference tests.
# Keep finite values on unrequired rows; only missing values are replaced by zero.
# Changing that convention can affect the canonical variance arithmetic.
.prepare_numeric <- function(time,delta,event,z,cc,time_transform=function(t)t) {
  is_integer <- function(x) is.numeric(x) && all(is.finite(x)) && all(x==trunc(x))
  if(!is.numeric(time) || !is.numeric(z) || !is.matrix(z))
    stop('X must be numeric and Z must be a numeric matrix.')
  for(x in list(delta,event,cc))
    if(!is_integer(x)) stop('delta, event_indic and cc must contain finite integers.')
  n <- length(time);p <- ncol(z)
  if(n==0L || nrow(z)!=n || length(delta)!=n || length(event)!=n || length(cc)!=n)
    stop('X, delta, event, cc and the rows of Z must have equal, non-zero lengths.')
  if(p==0L) stop('Z must have at least one covariate.')
  if(any(!is.finite(time)) || any(time<0)) stop('Follow-up times must be finite and non-negative.')
  if(any(!delta %in% 0:1) || any(!cc %in% 0:1) || any(event<0))
    stop('delta and cc must be 0/1; event must be a non-negative integer (0 means censored).')
  if(any((event==0)!=(delta==0))) stop('delta is inconsistent with event: event=0 must mean censored.')
  measured <- cc==1 | event==1
  if(any(!is.finite(z[measured,,drop=FALSE])))
    stop('Covariates are missing or non-finite for a subcohort member or target-event case.')
  if(any(is.infinite(z))) stop('Infinite covariates are not supported.')
  z[is.na(z)] <- 0
  if(!any(event==1)) stop('At least one target event (event=1) is required.')
  if(!any(cc==1)) stop('The subcohort must contain at least one subject.')
  centered <- scale(z[measured,,drop=FALSE],center=TRUE,scale=FALSE)
  sizes <- sqrt(colSums(centered^2))
  if(any(sizes==0)) stop('A covariate is constant among measured subjects; its effect is not identifiable.')
  centered <- sweep(centered,2,sizes,'/')
  if(qr(centered,tol=1e-10)$rank<p) stop('Covariates are linearly dependent among measured subjects.')
  if(!is.function(time_transform)) stop('time_transform must be a function.')
  q <- time_transform(time)
  if(is.numeric(q) && is.null(dim(q)) && length(q)==n) q <- matrix(q,n,p)
  if(!is.numeric(q) || !is.matrix(q) || !identical(dim(q),dim(z)) || any(!is.finite(q)))
    stop('time_transform must return finite numeric values, one per subject, or a matrix with the dimensions of Z.')
  if(is.null(colnames(z))) colnames(z) <- paste0('covariate',seq_len(p))
  list(time=as.numeric(time),delta=as.integer(delta),event=as.integer(event),z=z,
       cc=as.integer(cc),measured=measured,event_index=which(event==1),q=q,
       time_transform=time_transform)
}
