#' Title
#'
#' @param .x
#' @param pnum
#' @param panel
#'
#' @return
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
plot_donut <- function(patient, reference,  panel=NULL) {
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
  SummarizedExperiment::assays(patient)$abundance <- context_zero_markers(patient)

  sample <- SummarizedExperiment::assays(patient)$abundance |>
    tibble::as_tibble(rownames="Protein_Peptide") |>
    magrittr::set_colnames(c("Protein_Peptide","Patient"))


  ggplot2::ggplot(sample, ggplot2::aes(x=1, y=Patient, label=Protein_Peptide, fill=Protein_Peptide)) +
    ggplot2::geom_col() +
    ggplot2::coord_polar(theta="y") +
    #geom_text(
    #  aes(label = Protein_Peptide),
    #  position = position_stack(vjust = 0.5),
    #  size=8
    #) +
    ggplot2::theme_void() +
    ggplot2::scale_fill_brewer() +
    ggplot2::xlab("") + ggplot2::ylab("") +
    #theme(legend.title = element_blank()) +
    ggplot2::theme(legend.title=ggplot2::element_blank(),
                   plot.title=ggplot2::element_text(margin=ggplot2::margin(b=-50),vjust=1)) +

    #ggplot2::labs(title=main) +
    ggplot2::xlim(c(-0.5,1.5)) +
    # Put the text in the middle of the donut.
    ggplot2::annotate(geom = 'text', x = -0.5, y = 1, label = panel_name, size=4)

}

globalVariables(c("Patient","Protein_Peptide"))
