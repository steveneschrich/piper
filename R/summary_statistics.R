#' Title
#'
#' @param sample A sample to calculate relative to
#' @param reference A reference population
#' @param panels A list of panels to consider when calculating
#'
#' @return A table of summary statistics by assay
#' @export
#'
#' @importFrom methods is
#' @examples
summary_statistics_table <- function(sample, reference,  panels) {
    # Checking input parameters
    assertthat::assert_that(methods::is(reference, "SummarizedExperiment"))
    assertthat::assert_that(methods::is(sample, "SummarizedExperiment"))
    assertthat::assert_that(nrow(sample) == nrow(reference))

    # Panels needs to be one or more panels annotating the markers in reference/sample
    assertthat::assert_that(is.list(panels))
    assertthat::assert_that(length(panels)>0)


    # Panel: Select the specific markers to plot (all markers in sample)
    markers <- rownames(sample)

    # Subset the sample/reference data
    sample <- sample[markers,]
    reference <- reference[markers,]

    # Build the output table.
    output_table <- tibble::tibble(
      Markers = markers,
      # This appears to be make a list of panels that the marker is part of
      Panel = markers_in_panels(markers, panels),
      `Expression` = SummarizedExperiment::assay(sample)[,1],
      `Expression Percentile` = row_percent(sample, SummarizedExperiment::rowData(reference)$ecdf)[,1],
      `Reference Minimum Expression` = row_min(reference),
      `Reference Mean Expression` = row_mean(reference),
      `Reference Median Expression` = row_median(reference),
      `Reference Maximum Expression` = row_max(reference)
    )

  output_table
}





