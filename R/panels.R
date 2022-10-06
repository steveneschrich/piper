#' Title
#'
#' @param .x
#'
#' @return
#' @export
#'
#' @examples
panels <- function(.x) {
  unique(Biobase::fData(.x)$Subcategory)
}


#' Title
#'
#' @param .x
#' @param panel
#'
#' @return
#' @export
#' @importFrom rlang .data
#' @examples
get_panel <- function(.x, panel = "Tissue QC") {
  stopifnot(class(.x) %in% c("ExpressionSet"))

  # Extract specific sample (either the value or name, shouldn't matter)
  res <- Biobase::exprs(.x) |>
    data.frame() |>
    tibble::rownames_to_column("Protein_Peptide") |>
    # Add assay annotation
    dplyr::left_join(
      Biobase::fData(.x),
      by = "Protein_Peptide"
    ) |>
    dplyr::filter(.data$Subcategory %in% panel) |>
    dplyr::select(.data$Protein_Peptide, .data$`Protein Symbol`, .data$`Category`, .data$`Subcategory`, dplyr::everything()) |>
    tidyr::pivot_longer(
      cols = -tidyselect::any_of(c("Protein_Peptide","Protein Symbol","Category","Subcategory")),
      names_to="Patient",values_to="Expression")

  res
}

#' Title
#'
#' @param .x
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#' @examples
calculate_panel_statistics <- function(.x) {
  stopifnot(class(.x) %in% c("ExpressionSet"))

  # Compute max value per peptide
  Biobase::exprs(.x) |>
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
