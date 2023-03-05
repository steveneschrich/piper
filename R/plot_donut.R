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
plot_donut <- function(sample, reference,  panel=NULL) {
  # Checking input parameters
  assertthat::assert_that(is(reference, "SummarizedExperiment"))
  assertthat::assert_that(is(sample, "SummarizedExperiment"))
  assertthat::assert_that(nrow(sample) == nrow(reference))

  # Handle the panel input
  if (is.null(panel))
    panel <- list("Random"=base::sample(rownames(reference),5))
  assertthat::assert_that(is.list(panel))
  assertthat::assert_that(utils::hasName(panel,"markers"))


  # Panel: Select the specific panel to plot
  markers <- match_markers(panel$markers, reference)

  # Subset the sample/reference data
  sample <- sample[markers,]
  reference <- reference[markers,]

  # Statistics: We plot the percent of max abundance in the radar plot.
  SummarizedExperiment::assays(sample)$abundance <- impute_missing_with_values(sample, 0)


  sample <- SummarizedExperiment::assays(sample)$abundance |>
    tibble::as_tibble(rownames="Protein_Peptide") |>
    magrittr::set_colnames(c("Protein_Peptide","sample"))


  ggplot2::ggplot(sample, ggplot2::aes(x=1, y=sample, label=Protein_Peptide, fill=Protein_Peptide)) +
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
    ggplot2::annotate(geom = 'text', x = -0.5, y = 1, label = panel$name, size=4)

}

globalVariables(c("sample","Protein_Peptide"))
