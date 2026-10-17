#' @rdname cc_fit_methods
#' @param covariates Optional character vector of formula term or expanded
#'   coefficient names to plot, in the requested order. NULL plots all columns.
#'   Selecting a factor term includes all its expanded coefficients. Selection
#'   changes only the displayed plots, never the fitted model.
#' @export
plot.cc_fit <- function(x, covariates=NULL, ...) .plot_residuals(x,covariates)

.plot_residuals <- function(x,covariates=NULL) {
  labels <- colnames(x$event_residuals)
  selected <- seq_along(labels)
  if(!is.null(covariates)) {
    if(!is.character(covariates) || !length(covariates) || anyNA(covariates))
      stop("covariates must be a non-empty character vector of term or coefficient names.")
    terms <- x$coefficient_terms
    selected <- unique(unlist(lapply(covariates,function(name) {
      term_matches <- if(is.null(terms)) rep(FALSE,length(labels)) else terms==name
      matches <- which(labels==name | term_matches)
      if(!length(matches)) stop(paste("Unknown covariate or coefficient:",name))
      matches
    }),use.names=FALSE))
  }
  if (!isTRUE(x$converged)) {
    warning("Residual diagnostics are unavailable because the fit did not converge.")
    return(invisible(list()))
  }
  available <- vapply(selected,function(j) {
    r <- x$event_residuals[,j]
    length(r)>=2L && all(is.finite(r)) && all(is.finite(x$event_times)) &&
      length(unique(r))>1L && length(unique(x$event_times))>1L
  },logical(1))
  if(any(!available))
    warning(paste('Residual plot unavailable for',paste(labels[selected[!available]],collapse=', '),
                  ': non-finite values or insufficient variation.'),call.=FALSE)
  selected <- selected[available]
  if(!length(selected))return(invisible(list()))
  Time <- Residuals <- NULL # Names evaluated inside ggplot2::aes.
  order <- order(x$event_times, x$event_index)
  plots <- lapply(selected, function(j) {
    df <- data.frame(Time = x$event_times[order], Residuals = x$event_residuals[order, j])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Residuals)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
      ggplot2::geom_point(size=0.9,alpha=0.5) +
      ggplot2::labs(title = colnames(x$event_residuals)[j], x = "Time", y = "Schoenfeld-type residual") +
      ggplot2::theme_minimal()
    if (nrow(df) >= 10L && length(unique(df$Time)) >= 5L)
      p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "red")
    print(p)
    p
  })
  names(plots) <- labels[selected]
  invisible(plots)
}
