
#' Title
#'
#' @param .x
#' @param pnum
#' @param panel
#' @param use_percent
#' @param peptide_list
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @examples
plot_radar <- function(.x, pnum = 1,
                        panel = "Tissue QC",
                        use_percent = FALSE,
                        peptide_list=NULL,
                        axis.label.size=5) {

  stopifnot(class(.x) %in% c("ExpressionSet"))

  if ( is.character(pnum))
    pnum <- which(Biobase::sampleNames(.x) %in% pnum)

  qc <- get_sample_statistics(.x, pnum, panel = panel) %>%
    # Reduce the total to only the two needed fields (Peptide and Percent)
    dplyr::select(.data$Protein_Peptide, .data$Percent)

  # Filter on the ones requested by peptide_list, if provided
  if (!is.null(peptide_list)) {
    qc <- qc %>% dplyr::filter(.data$Protein_Peptide %in% peptide_list)
  }
  qc <- qc %>%
    # Clean assay name
    dplyr::mutate(Protein_Peptide = clean_assay_name(.data$Protein_Peptide)) %>%
    # Transpose for wider format.
    tidyr::pivot_wider(names_from=.data$Protein_Peptide, values_from=.data$Percent) %>%
    # radar plot needs a first column that is not data
    dplyr::mutate(Patient = pnum) %>%
    dplyr::select(.data$Patient, dplyr::everything())

  ggradar::ggradar(qc, axis.label.size=axis.label.size) +
    ggplot2::ggtitle(panel)

}

# NB: res$x[[1]]$hovertext
