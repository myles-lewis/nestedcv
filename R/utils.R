
#' @importFrom utils argsAnywhere

checkArg <- function(x, fun) {
  args <- formalArgs(fun)
  if (any(args == x)) return(TRUE)
  # for S3/4 methods
  met <- suppressWarnings(methods(fun))
  args <- unlist(lapply(met, function(i) formalArgs(argsAnywhere(i))))
  any(args == x)
}


#' Supervised PCA plot
#' 
#' Performs supervised principle component analysis (PCA) after filtering
#' dataset to help determine whether filtering has been useful for separating
#' samples according to the outcome variable.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter].
#'   Any function can be provided and is passed `y` and `x`. Must return a
#'   character vector with names of filtered predictors.
#' @param filter_options List of additional arguments passed to the filter
#'   function specified by `filterFUN`.
#' @param plot Logical whether to plot a ggplot2 object or return the PC scores
#' @param ... Optional arguments passed to [princomp()]
#' @return If `plot=TRUE` returns a ggplot2 plot, otherwise returns the
#'   principle component scores.
#' 
#' @importFrom stats princomp
#' @importFrom ggplot2 theme_classic
#' @export

supervisedPCA <- function(y, x,
                          filterFUN = NULL,
                          filter_options = NULL,
                          plot = TRUE, ...) {
  ok <- checkxy(y, x)
  y <- y[ok$r]
  x <- x[ok$r, ok$c]
  dat <- nest_filt_bal(NULL, y, x, filterFUN, filter_options)
  filtx <- dat$filt_xtrain
  df <- princomp(filtx, ...)
  data <- data.frame(df$scores)
  data$outcome <- y
  if (!plot) return(data)
  ggplot(data, aes(x = .data$Comp.1, y = .data$Comp.2, color = .data$outcome)) +
    geom_point() +
    theme_classic()
}


inspectGrid <- function(method, tuneLength = 3, y = NA, x = NA) {
  tgrid <- caret::getModelInfo(method)[[1]]$grid(x, y, tuneLength)
  cat("Tuning grid with", ncol(tgrid), "parameter(s) and", nrow(tgrid), "rows\n")
  choices <- lapply(tgrid, unique)
  maxchar <- max(unlist(lapply(names(choices), nchar)))
  lapply(names(choices), function(param) {
    cat(param, spaces(maxchar - nchar(param) +2))
    cat(paste(choices[[param]], collapse=", "), "\n")
  })
  invisible(choices)
}

spaces <- function(n) {
  paste0(rep(" ", n), collapse = "")
}
