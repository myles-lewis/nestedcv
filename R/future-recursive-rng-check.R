#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom stats runif
check_future_recursive_rng <- function() {
    oplan <- plan(future::sequential)
    on.exit(plan(oplan))
    y <- do.call(rbind, future_lapply(1:3, future.seed = 1986, FUN = function(i) {
        do.call(rbind, future_lapply(1:3, future.seed = TRUE, FUN = function(j) {
            data.frame(i = i, j = j, random = runif(n = 1L)) }))
    }))
    if (anyDuplicated(y$random)) {
        warning(
          "Your installed copy of future.apply has a bug that causes RNG to produce repeated values in nested future_lapply calls (see future.apply issue #108). This package is likely to produce invalid and misleading results. A version of the future.apply with a fix merged is available here: https://github.com/DarwinAwardWinner/future.apply. You can install it with devtools::install_github('DarwinAwardWinner/future.apply').",
          call. = FALSE
        )
    }
    invisible()
}

.onLoad <- function(libname, pkgname) {
    check_future_recursive_rng()
}
