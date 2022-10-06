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
#'
#' @examples
#' \dontrun{
#' plot_rug(lscc, pnum = 1, panel = "Tissue QC")
#' }
plot_rug <- function(.x, pnum=1, panel="Tissue QC") {
  if ( is.character(pnum))
    pnum <- which(Biobase::sampleNames(.x) %in% pnum)

  ip <- which(Biobase::fData(.x)$Subcategory == panel)

  min_peptide_level <- apply(Biobase::exprs(.x), 1, min,na.rm=T)-0.1

  sample <- Biobase::exprs(.x)[,pnum]
  sample <- ifelse(is.na(sample), min_peptide_level, sample) |>
    tibble::enframe(name="Protein_Peptide", value="Patient") |>
    dplyr::left_join(Biobase::fData(.x),
                     by=c("Protein_Peptide"="Protein_Peptide"))

  pl <- sapply(ip, function(n) {

    mydf <- data.frame(x=tidyr::replace_na(Biobase::exprs(.x)[n,]),0)
    p <- ggplot2::ggplot(mydf,ggplot2::aes(x=x)) +
      ggplot2::geom_rug(length=grid::unit(1,"npc")) +
      ggplot2::coord_cartesian(ylim=c(0,1),clip = 'off') +

      ggplot2::geom_point(ggplot2::aes(x=sample$Patient[n],y=c(-0.01)),shape=24, fill="red",size=4) +
      ggplot2::geom_point(ggplot2::aes(x=sample$Patient[n],y=c(1.01)), shape=25,
                          fill="red",size=4) +
      ggplot2::xlab("") +
      ggplot2::ylab(Biobase::featureNames(.x)[n]) +
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
      title = panel,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face="bold"),
                             plot.background=ggplot2::element_rect(size=2, color="black")))

  list(plot_panel)

}

globalVariables("x")
