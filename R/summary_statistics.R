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
summary_statistics_all <- function(patient, reference,  panels) {
    # Checking input parameters
    assertthat::assert_that(is(reference, "SummarizedExperiment"))
    assertthat::assert_that(is(patient, "SummarizedExperiment"))
    assertthat::assert_that(nrow(patient) == nrow(reference))

    # Panels needs to be one or more panels annotating the markers in reference/patient
    assertthat::assert_that(is.list(panels))
    assertthat::assert_that(length(panels)>0)

    # Create a panel annotation file
    panel_annotation <- panels |>
      tibble::enframe(name = "Panel", value = "Marker") |>
      tidyr::unchop(Marker)


    # Panel: Select the specific markers to plot (all markers in sample)
    markers <- rownames(patient)

    # Subset the patient/reference data
    patient <- patient[markers,]
    reference <- reference[markers,]

    # Build the output table.
    output_table <- tibble::tibble(
      Markers = markers,
      Panel = annotate_markers_with_panels(markers, panel_annotation),
      `Expression` = SummarizedExperiment::assay(patient)[markers,],
      `Expression Quantile` = purrr::map_dbl(markers, function(m) {
        stats::ecdf(SummarizedExperiment::assay(reference)[m,])(
          SummarizedExperiment::assay(patient)[m,]
        )
      }),
      t(
        apply(SummarizedExperiment::assay(reference), 1, function(v) {
          c(
            `Reference Minimum Expression` = min(v,na.rm=T),
            `Reference Mean Expression` = mean(v, na.rm=T),
            `Reference Maximum Expression` = max(v, na.rm=T)
          )
        })
      )
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
