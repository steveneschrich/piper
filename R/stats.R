# context.R
#' Calculate percent of values for samples
#'
#' @param x A [SummarizedExperiment::SummarizedExperiment()] representing samples.
#' @param values Values to use as the maximum/total from which to calculate percentages.
#' @param use_na If TRUE (default), keep NA values from sample as NA percent. Otherwise, use 0 percent.
#'
#' @return A fractional percentage calculated against `values` for the default `assay` slot
#'  of `x`.
#' @export
#'
#' @examples
#' \dontrun{
#' irisse <- SummarizedExperiment::SummarizedExperiment(t(iris[,1:4]))
#' context_percent_max(irisse[,5], irisse)
#' }
old_percent <- function(x, values, use_na = TRUE) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))
  assertthat::assert_that(nrow(x) == length(values))

  # Percentage of value is a simple calculation here.
  p <- SummarizedExperiment::assay(x) / values

  if ( !use_na )
    p <- tidyr::replace_na(p, 0)

  p
}

#' Title
#'
#' @param x
#' @param ecdfs
#'
#' @return
#' @export
#'
#' @examples
row_percent <- function(x, ecdfs, use_na = TRUE) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))
  assertthat::assert_that(is.list(ecdfs))
  assertthat::assert_that(nrow(x)==length(ecdfs) && all(rownames(x)==names(ecdfs)))

  res <- apply(SummarizedExperiment::assay(x), 2, function(v) {
    purrr::map2_dbl(v, ecdfs, ~.y(.x))
  })
  if (!use_na)
    res <- apply(res, 2, function(.x) {tidyr::replace_na(.x, replace = 0)})


  res
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



#' Impute missing values with specified alternate values
#'
#' @param x A [SummarizedExperiment::SummarizedExperiment()] to impute missing values from.
#' @param values A vector of values (recycle rules per [dplyr::coalesce()]) to use as imputed values.
#' @param offset A value to add to values (default 0).
#'
#' @return A matrix the same dimensions as the default `assay` slot
#' changed such that missing values are imputed with the specified alternate values.
#'
#' @export
#'
#' @examples
impute_missing_with_values <- function(x, values, offset = 0) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))

  values <- values + offset

  xp <- apply(SummarizedExperiment::assay(x), 2, function(.x) {
    dplyr::coalesce(.x, values)
  })

  xp
}

#' Calculate row-wise minimums
#'
#' @param x A [SummarizedExperiment::SummarizedExperiment()] to calculate row minimums
#' @param na.rm Logical (default TRUE). Should NA values be removed from consideration (TRUE) or
#' treated as a 0 (FALSE).
#'
#' @return A named vector with each element the minimum of the corresponding `x` row. Each
#' element is named with the corresponding row name of `x`.
#'
#' @export
#'
#' @examples
row_min <- function(x, na.rm = TRUE) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))

  xv <- SummarizedExperiment::assay(x)

  if ( !na.rm )
    xv[is.na(xv)] <- 0

  matrixStats::rowMins(xv, na.rm = TRUE, useNames = TRUE)

}


#' Calculate row-wise maximums
#'
#' @param x A [SummarizedExperiment::SummarizedExperiment()] to calculate row maxmimums
#' @param na.rm Logical (default TRUE). Should NA values be removed from consideration (TRUE) or
#' treated as maxval (FALSE).
#' @param maxval A maximum value to use if there are NA's. Default is 0.
#'
#' @return A named vector with each element the maximum of the corresponding `x` row. Each
#' element is named with the corresponding row name of `x`.
#'
#' @export
#'
#' @examples
row_max <- function(x, na.rm = TRUE, maxval = 0) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))

  xv <- SummarizedExperiment::assay(x)
  if ( !na.rm )
    xv[is.na(xv)] <- maxval

  matrixStats::rowMaxs(xv, useNames = TRUE, na.rm = TRUE)

}

#' Title
#'
#' @param x
#' @param na.rm
#'
#' @return
#' @export
#'
#' @examples
row_mean <- function(x, na.rm = TRUE) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))

  xv <- SummarizedExperiment::assay(x)

 rowMeans(xv, na.rm = na.rm)

}

#' Title
#'
#' @param x
#' @param na.rm
#'
#' @return
#' @export
#'
#' @examples
row_median <- function(x, na.rm = TRUE) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))

  xv <- SummarizedExperiment::assay(x)

  matrixStats::rowMedians(xv, useNames = TRUE, na.rm = na.rm)

}


#' Calculate row-wise empirical cumulative distribution functions
#'
#' @param x A [SummarizedExperiment::SummarizedExperiment()] to calculate ecdfs from.
#'
#' @return A list of ecdfs
#' @export
#'
#' @examples
row_ecdf <- function(x) {
  assertthat::assert_that(methods::is(x, "SummarizedExperiment"))

  apply(SummarizedExperiment::assay(x), 1, function(v) {
    ecdf(v)
  })

}
