


#' @export
#' 
boot_filter <- function(y, x, filterFUN, B = 50, type = "index", ...) {
  ranks <- sapply(1:B, function(i) {
    ind <- sample.int(length(y), replace = TRUE)
    out <- filterFUN(y[ind], x[ind, ], type = "full", ...)
    rank(out[, grep("p.?val", colnames(out))])
  })
  meanRank <- rowMeans(ranks)
  names(meanRank) <- colnames(x)
  if (type == "full") return(sort(meanRank))
  order(meanRank)
}


#' @export
boot_ttest <- function(y, x, B = 50, ...) {
  boot_filter(y, x, ttest_filter, B=B, ...)
}

