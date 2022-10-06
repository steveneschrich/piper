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
plot_donut <- function(.x, pnum = 1, panel = "Tissue QC") {

  if ( is.character(pnum))
    pnum <- which(Biobase::sampleNames(.x) %in% pnum)


  qc <- get_sample(.x, pnum, use_na = FALSE) %>%
    # Extract out the panel of interest.
    dplyr::filter(.data$Subcategory %in% panel) %>%
    # Clean assay name
    dplyr::mutate(Protein_Peptide = clean_assay_name(.data$Protein_Peptide))

  ggplot2::ggplot(qc, ggplot2::aes(x=1, y=Patient, label=Protein_Peptide, fill=Protein_Peptide)) +
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
    ggplot2::annotate(geom = 'text', x = -0.5, y = 1, label = panel, size=4)

}

globalVariables(c("Patient","Protein_Peptide"))
