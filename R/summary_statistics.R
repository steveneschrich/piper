#' Title
#'
#' @param .x
#' @param pnum
#' @param panel
#'
#' @return
#' @export
#'
#' @examples
summary_statistics_table <- function(sample, reference,  panels) {
    # Checking input parameters
    assertthat::assert_that(is(reference, "SummarizedExperiment"))
    assertthat::assert_that(is(sample, "SummarizedExperiment"))
    assertthat::assert_that(nrow(sample) == nrow(reference))

    # Panels needs to be one or more panels annotating the markers in reference/sample
    assertthat::assert_that(is.list(panels))
    assertthat::assert_that(length(panels)>0)

    # Create a panel annotation file
    panel_annotation <- panels |>
      tibble::enframe(name = "Panel", value = "Marker") |>
      tidyr::unchop(Marker)


    # Panel: Select the specific markers to plot (all markers in sample)
    markers <- rownames(sample)

    # Subset the sample/reference data
    sample <- sample[markers,]
    reference <- reference[markers,]

    # Build the output table.
    output_table <- tibble::tibble(
      Markers = markers,
      Panel = annotate_markers_with_panels(markers, panel_annotation),
      `Expression` = SummarizedExperiment::assay(sample),
      `Expression Percentile` = row_percent(sample, SummarizedExperiment::rowData(reference)$ecdf),
      `Reference Minimum Expression` = row_min(reference),
      `Reference Mean Expression` = row_mean(reference),
      `Reference Median Expression` = row_median(reference),
      `Reference Maximum Expression` = row_max(reference)
    )

  output_table
}

#' Title
#'
#' @param markers
#' @param panels
#'
#' @return
#' @export
#'
#' @examples
annotate_markers_with_panels <- function(markers, panels) {
  assertthat::assert_that(is.character(markers))
  assertthat::assert_that(is.data.frame(panels))
  assertthat::assert_that(hasName(panels, c("Panel","Marker")))

  dplyr::left_join(
    tibble::enframe(markers, name=NULL, value = "Marker"),
    panels,
    by = "Marker"
  ) |>
    tidyr::chop(Panel) |>
    dplyr::pull(Panel) |>
    purrr::map_chr(~stringr::str_c(.x, collapse=", "))
}
