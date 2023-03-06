#' Evaluate a score for a panel
#'
#' @param sample
#' @param reference
#' @param panel
#'
#' @return
#' @export
#'
#' @examples
eval_score <- function(sample, reference, panel) {
  if (!utils::hasName(panel, "scorefn"))
    sf <- "piper::max_percentile"
  else
    sf <- panel$scorefn

  pkg_fn <- stringi::stri_split_fixed(sf, "::",2)[[1]]

  fn <- get(pkg_fn[2], envir=asNamespace(pkg_fn[1]))
  fn(sample, reference, panel)
}


#' Title
#'
#' @param sample
#' @param reference
#' @param panel
#'
#' @return
#' @export
#'
#' @examples
max_percentile <- function(sample, reference, panel) {
  percentile(sample, reference, panel) |>
    max()
}

#' Title
#'
#' @param sample
#' @param reference
#' @param panel
#'
#' @return
#' @export
#'
#' @examples
median_percentile <- function(sample, reference, panel) {
  percentile(sample, reference, panel) |>
    median()
}
#' Title
#'
#' @param sample
#' @param reference
#' @param panel
#'
#' @return
#' @export
#'
#' @examples
min_percentile <- function(sample, reference, panel) {
  percentile(sample,reference, panel) |>
    min()
}

#' Title
#'
#' @param sample
#'
#' @return
#' @export
#'
#' @examples
percent_missing <- function(sample, reference, panel) {
  nm <- length(
    which(
      !is.na(
        SummarizedExperiment::assay(sample)[,1]
      )
    )
  )
  nm/nrow(sample)
}
