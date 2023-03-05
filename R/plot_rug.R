#' Generate a rug plot from proteomics panel
#'
#' @description Provide slider-style readout of proteomics panel
#'  components.
#'
#' @param sample A SummarizedExperiment representing a single sample
#' @param reference A SummarizedExperiment representing the reference population
#' @param panel The panel to visualize (list consisting of name and markers).
#'
#' @return A [patchwork::patchwork()] panel of graphs representing a rug plot.
#' @export
#' @examples
#' \dontrun{
#' plot_rug(lscc[,1], lscc , panel = panels$ADC)
#' }
plot_rug <- function(sample, reference,  panel=NULL) {
  # Checking input parameters
  assertthat::assert_that(methods::is(reference, "SummarizedExperiment"))
  assertthat::assert_that(methods::is(sample, "SummarizedExperiment"))
  assertthat::assert_that(nrow(sample) == nrow(reference))

  assertthat::assert_that(is.list(panel))
  assertthat::assert_that(utils::hasName(panel, "markers"))

  # Panel: Select the specific panel to plot
  markers <- match_markers(panel$markers, reference)
  panel_name <- panel$name

  # Subset the sample/reference data
  sample <- sample[markers,]
  reference <- reference[markers,]

  # Statistics: We plot the percent of max abundance in the radar plot. We use the minimum value of the reference
  # set (offset by -0.1, to show it is lower than anything else) if the sample value is missing.
  SummarizedExperiment::assays(sample)$abundance <- impute_missing_with_values(
    sample,
    values = row_min(reference),
    offset = -0.1
  )



  sample_table <- SummarizedExperiment::assays(sample)$abundance |>
    tibble::as_tibble(rownames="Protein_Peptide") |>
    magrittr::set_colnames(c("Protein_Peptide","sample")) |>
    tibble::deframe()

  pl <- sapply(names(sample_table), function(n) {

    mydf <- data.frame(x=tidyr::replace_na(SummarizedExperiment::assay(reference)[n,],0))
    p <- ggplot2::ggplot(mydf,ggplot2::aes(x=x)) +
      ggplot2::geom_rug(length=grid::unit(1,"npc")) +
      ggplot2::coord_cartesian(ylim=c(0,1),clip = 'off') +

      ggplot2::geom_point(ggplot2::aes(x=sample_table[n],y=c(-0.01)),shape=24, fill="red",size=4) +
      ggplot2::geom_point(ggplot2::aes(x=sample_table[n],y=c(1.01)), shape=25,
                          fill="red",size=4) +
      ggplot2::xlab("") +
      ggplot2::ylab(n) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.background=ggplot2::element_rect(),
                     panel.grid = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust=1)
      )
    list(p)
  })

  plot_panel <- patchwork::wrap_plots(pl,ncol=1) +
    patchwork::plot_annotation(
      title = panel_name,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face="bold"),
                             plot.background=ggplot2::element_rect(size=2, color="black")))

  list(plot_panel)

}

globalVariables("x")
