#' Generate a rug plot from proteomics panel
#'
#' @description Provide slider-style readout of proteomics panel
#'  components.
#'
#' @param .x The [Biobase::ExpressionSet] to visualization
#' @param pnum For multi-patient datasets, the patient to visualize
#' (should be in `colnames(.x)` or a numeric index/column number).
#' @param panel The panel to visualize
#'
#' @return A [patchwork::patchwork()] panel of graphs representing a rug plot.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' plot_rug(lscc, pnum = 1, panel = "Tissue QC")
#' }
plot_rug <- function(patient, reference,  panel=NULL) {
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
  SummarizedExperiment::assays(patient)$abundance <- context_floored_markers(patient, reference)


  sample <- SummarizedExperiment::assays(patient)$abundance |>
    tibble::as_tibble(rownames="Protein_Peptide") |>
    magrittr::set_colnames(c("Protein_Peptide","Patient")) |>
    tibble::deframe()

  pl <- sapply(names(sample), function(n) {

    mydf <- data.frame(x=tidyr::replace_na(SummarizedExperiment::assay(reference)[n,],0))
    p <- ggplot2::ggplot(mydf,ggplot2::aes(x=x)) +
      ggplot2::geom_rug(length=grid::unit(1,"npc")) +
      ggplot2::coord_cartesian(ylim=c(0,1),clip = 'off') +

      ggplot2::geom_point(ggplot2::aes(x=sample[n],y=c(-0.01)),shape=24, fill="red",size=4) +
      ggplot2::geom_point(ggplot2::aes(x=sample[n],y=c(1.01)), shape=25,
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
