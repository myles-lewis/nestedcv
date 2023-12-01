#' Text Progress Bar 2
#' 
#' Text progress bar in the R console. Modified from `utils::txtProgressBar()`
#' to include title and timing.
#' 
#' @param min Numeric value for minimum of the progress bar.
#' @param max Numeric value for maximum of the progress bar.
#' @param initial Initial value for the progress bar.
#' @param char The character (or character string) to form the progress bar.
#' @param width The width of the progress bar, as a multiple of the width of
#'   `char`. If `NA`, the default, the number of characters is that which fits
#'   into `getOption("width")`.
#' @param title Title for the progress bar.
#' @details
#' Use `utils::setTxtProgressBar()` to set the progress bar and `close()` to
#' close it.
#' 
#' @importFrom utils flush.console
#' @returns An object of class "`txtProgressBar`".
#' @export

txtProgressBar2 <- function(min = 0, max = 1, initial = 0, char = "=",
                            width = NA, title = "") {
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L
  nw <- nchar(char, "w")
  if (nw == 0) 
    stop("'char' must have a non-zero width")
  if (is.na(width)) {
    width <- getOption("width")
    nt <- nchar(title)
    if (length(nt) == 0) nt <- 0
    width <- width - 22L - nt
    if (nw > 1) 
      width <- trunc(width/nw)
  }
  if (max <= min) 
    stop("must have 'max' > 'min'")
  up3 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    cat(paste0("\r", title, "  |", strrep(" ", nw * width + 6)))
    cat(paste(c("\r", title, "  |", rep.int(char, nb),
                rep.int(" ", nw * (width - nb)),
                sprintf("| %3d%%", pc)), collapse = ""))
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() if (!.killed) {
    end <- Sys.time()
    cat(paste0("  (", format(end - .start, digits = 3), ")\n"))
    flush.console()
    .killed <<- TRUE
  }
  up <- up3
  up(initial)
  .start <- Sys.time()
  structure(list(getVal = getVal, up = up, kill = kill),
            class = "txtProgressBar")
}
