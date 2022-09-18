
#' @importFrom utils argsAnywhere

checkArg <- function(x, fun) {
  args <- formalArgs(fun)
  if (any(args == x)) return(TRUE)
  # for S3/4 methods
  met <- suppressWarnings(methods(fun))
  args <- unlist(lapply(met, function(i) formalArgs(argsAnywhere(i))))
  any(args == x)
}

