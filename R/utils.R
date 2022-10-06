
#' @importFrom utils argsAnywhere

checkArg <- function(x, fun) {
  args <- formalArgs(fun)
  if (any(args == x)) return(TRUE)
  # for S3/4 methods
  met <- suppressWarnings(methods(fun))
  args <- unlist(lapply(met, function(i) formalArgs(argsAnywhere(i))))
  any(args == x)
}



#' @export

supervisedPCA <- function(y, x,
                          filterFUN = NULL,
                          filter_options = NULL,
                          plot = TRUE, ...) {
  dat <- nest_filt_bal(NULL, y, x, filterFUN, filter_options)
  filtx <- dat$filt_xtrain
  df <- princomp(filtx, ...)
  data <- data.frame(df$scores)
  data$outcome <- y
  if (!plot) return(data)
  ggplot(data, aes(x = Comp.1, y = Comp.2, color = outcome)) +
    geom_point() +
    theme_classic()
}

