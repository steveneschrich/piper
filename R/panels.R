#' Title
#'
#' @param .x
#'
#' @return
#' @export
#'
#' @examples
panels <- function(.x) {
  unique(SummarizedExperiment::rowData(.x)$Subcategory)
}


#' Title
#'
#' @param .x
#' @param panel
#'
#' @return
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @examples
get_panel <- function(.x, panel = "Tissue QC") {
  stopifnot(methods::is(.x, "SummarizedExperiment"))

  # Extract specific sample (either the value or name, shouldn't matter)
  res <- SummarizedExperiment::assay(.x) |>
    data.frame() |>
    tibble::rownames_to_column("Protein_Peptide") |>
    # Add assay annotation
    dplyr::left_join(
      SummarizedExperiment::rowData(.x),
      by = "Protein_Peptide"
    ) |>
    dplyr::filter(.data$Subcategory %in% panel) |>
    dplyr::select(.data$Protein_Peptide, .data$`Protein Symbol`, .data$`Category`, .data$`Subcategory`, dplyr::everything()) |>
    tidyr::pivot_longer(
      cols = -tidyselect::any_of(c("Protein_Peptide","Protein Symbol","Category","Subcategory")),
      names_to="Patient",values_to="Expression")

  res
}

#' Calculate statistics for a given panel
#'
#' @param .x
#' @param panel A assay panel, or a list with (name,markers)
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#' @examples
calculate_panel_statistics <- function(.x, panel) {
  stopifnot(methods::is(.x, "SummarizedExperiment"))
  stopifnot(utils::hasName(panel, "markers"))

  # Compute max value per peptide
  SummarizedExperiment::assay(.x) |>
    # exprs is a matrix, so convert to data frame.
    data.frame() |>
    # move rownames to a new column
    tibble::rownames_to_column("Protein_Peptide") |>
    # All observations are an individual row
    tidyr::pivot_longer(
      cols=c(dplyr::everything(), -.data$Protein_Peptide),
      names_to="Sample",
      values_to="Expression"
    ) |>
    # Group by the peptide
    dplyr::group_by(.data$Protein_Peptide) |>
    # Lots of summaries
    dplyr::summarize(
      MaxExpression = max(.data$Expression, na.rm = TRUE),
      MinExpression = min(.data$Expression, na.rm = TRUE),
      MeanExpression = mean(.data$Expression, na.rm = TRUE),
      MedianExpression = stats::median(.data$Expression, na.rm = TRUE),
      ECDF = list(stats::ecdf(.data$Expression)),
      .groups="drop"
    )


}




#' Determine if marker is in a list of panels
#'
#' @description A list of panels each contains a set of markers. Test
#' if a given marker is in each of the panels. Returns a (named) logical vector
#' indicating if the marker is the panel.
#'
#' @param marker A string representing a marker
#' @param panels A list of panels, each of which has a markers component.
#'
#' @return A list of named logicals, with the name representing the
#'  panel and the value representing if the marker is within the panel.
#' @export
#'
#' @examples
#' \dontrun{
#' markers_in_panels("CDN2A_ALLEAGALPNAPNSYGR", panels)
#' }
is_marker_in_panels <- function(marker, panels) {
  purrr::map_lgl(panels, \(p) {
    marker %in% p[["markers"]]
  })
}

#' Title
#'
#' @param marker
#' @param panels
#'
#' @return
#' @export
#'
#' @examples
marker_in_panels <- function(marker, panels) {
  in_panels <- is_marker_in_panels(marker, panels)
  names(panels)[in_panels] |>
    paste(collapse = ", ")
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
markers_in_panels <- function(markers, panels) {
  assertthat::assert_that(is.character(markers))
  assertthat::assert_that(is.list(panels))

  purrr::map_chr(markers, \(m) marker_in_panels(m, panels))
}

