
#' Title
#'
#' @param .x
#' @param patient
#' @param panel
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @examples
plot_radar <- function(patient, reference, panel=NULL,
                       axis.label.size=4) {

  # Checking input parameters
  assertthat::assert_that(is(reference, "SummarizedExperiment"))
  assertthat::assert_that(is(patient, "SummarizedExperiment"))
  assertthat::assert_that(nrow(patient) == nrow(reference))

  # Handle the panel input
  if (is.null(panel))
    panel <- list("Random"=sample(rownames(reference),5))
  assertthat::assert_that(is.list(panel))
  assertthat::assert_that(length(panel[[1]])>0)
  assertthat::assert_that(all(is.numeric(panel[[1]])) || all(panel[[1]] %in% rownames(reference)))

  # Panel: Select the specific panel to plot
  markers <- match_markers(panel[[1]], reference)
  panel_name <- names(panel)[1]

  # Subset the patient/reference data
  patient <- patient[markers,]
  reference <- reference[markers,]

  # Statistics: We plot the percent of max abundance in the radar plot.
  SummarizedExperiment::assays(patient)$percent_max <- context_percent_max(patient, reference, use_na = FALSE)


  # Transform to plot data.
  qc <- SummarizedExperiment::assays(patient)$percent_max |>
    tibble::as_tibble(rownames = "Protein_Peptide") |>
    magrittr::set_colnames(c("Protein_Peptide","Percent")) |>
    # Transpose to wider format
    tidyr::pivot_wider(names_from=Protein_Peptide, values_from=Percent) |>
    # radar plot needs a first column that is not data
    dplyr::mutate(Patient = dplyr::row_number()) |>
    dplyr::select(Patient, dplyr::everything())

  ggradar::ggradar(qc, axis.label.size=axis.label.size) +
    ggplot2::ggtitle(panel_name)

}

# NB: res$x[[1]]$hovertext


