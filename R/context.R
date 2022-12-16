# context.R
#' Sample abundance as percentage of max reference abundance
#'
#' @param sample A [SummarizedExperiment::SummarizedExperiment()] representing a single sample
#' @param reference A [SummarizedExperiment::SummarizedExperiment()] representing the entire sample population
#' @param use_na If TRUE (default), keep NA values from sample as NA percent. Otherwise, use 0 percent.
#'
#' @return A one-column matrix of percent max (fractional)
#' @export
#'
#' @examples
#' \dontrun{
#' irisse <- SummarizedExperiment::SummarizedExperiment(t(iris[,1:4]))
#' context_percent_max(irisse[,5], irisse)
#' }
context_percent_max <- function(sample, reference, use_na = TRUE) {
  assertthat::assert_that(
    is_comparable(sample,reference),
    msg = "Sample and reference are not comparable, cannot calculate percent max."
  )

  #  check_stuff(patient, rerference)
  maxvals <- SummarizedExperiment::assay(sample) / max(SummarizedExperiment::assay(reference), na.rm=TRUE)

  if ( !use_na )
    maxvals <- tidyr::replace_na(maxvals, 0)

  maxvals
}

#' Title
#'
#' @param reference
#' @param sample
#'
#' @return
#' @export
#'
#' @examples
context_floored_markers <- function(sample, reference) {
  assertthat::assert_that(ncol(sample)==1)

  min_levels <- apply(SummarizedExperiment::assay(reference), 1, min, na.rm = TRUE) - 0.1
  apply(SummarizedExperiment::assay(sample), 2, function(.x) {
    dplyr::coalesce(.x, min_levels)
  })
}

#' Title
#'
#' @param sample
#'
#' @return
#' @export
#'
#' @examples
context_zero_markers <- function(sample) {
  assertthat::assert_that(ncol(sample)==1)

  apply(SummarizedExperiment::assay(sample), 2, function(.x) {
    dplyr::coalesce(.x, 0)
  })
}

#' Test if sample and reference are comparable
#'
#' @param sample
#' @param reference A [SummarizedExperiment::SummarizedExperiment()] representing the reference population.
#'
#' @return Logical value indicating if sample and reference are comparable.
#' @export
#'
#' @examples
is_comparable <- function(sample, reference) {

  overall_status <- TRUE


  if ( !nrow(sample) == nrow(reference) )  {
    logger::log_info("Number of measurements in sample (p={nrow(sample)}) not equal to reference population (p={nrow(reference)}).")
    overall_status <- FALSE
  }

  if (! all(rownames(sample) %in% rownames(reference))) {
    logger::log_info("Assay names are not consistent between sample and reference, cannot be merged.")
    overall_status <- FALSE
  }

  if (! all(rownames(sample) == rownames(reference))) {
    logger::log_info("Assay names are not matched (in order) between sample and reference, cannot be merged.")
    overall_status <- FALSE
  }

  overall_status
}


min_expression <- function(x, na_as_zero = FALSE) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))

  apply(SummarizedExperiment::assay(x), 1, function(v) {
    min(v, na.rm = TRUE)
  })
}
