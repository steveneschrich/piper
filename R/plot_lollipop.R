#' Create a lollipop plot of markers vs reference
#'
#' @param sample A sample SummarizedExperiment
#' @param reference A reference population
#' @param panel A panel (list with name and markers)
#'
#' @return
#' @export
#' @examples
plot_lollipop <- function(sample, reference, panel) {

  # Checking input parameters
  assertthat::assert_that(is(reference, "SummarizedExperiment"))
  assertthat::assert_that(is(sample, "SummarizedExperiment"))
  assertthat::assert_that(nrow(sample) == nrow(reference))

  assertthat::assert_that(is.list(panel))
  assertthat::assert_that(utils::hasName(panel,"markers"))

  # Panel: Select the specific panel to plot
  markers <- match_markers(panel$markers, reference)

  # Subset the sample/reference data
  sample <- sample[markers,]
  reference <- reference[markers,]

  patient_name <- colnames(sample)[1]

  rev_label <- function(s) {
    x <- stringr::str_split_fixed(s, "_",n= 2)
    paste0(x[,2],"_",x[,1])
  }


    # Note sure why this is here
    #dplyr::filter(Subcategory %in% c("Cancer Antigens","CDK","MAPK","RTK","PI3K/AKT/MTOR","Tissue QC")) |>
    dplyr::mutate(Protein_Peptide=stringr::str_replace_all(.data$Protein_Peptide, "_","\n"))

    gdf <- tibble::tibble(
      assay_id = rownames(reference),
      assay_label = stringr::str_replace(rev_label(.data$assay_id), "_","\n"),
      min = row_min(reference),
      max = row_max(reference),
      target = tidyr::replace_na(SummarizedExperiment::assay(sample)[,1],0),
      percentile =  row_percent(
        sample,
        ecdfs = SummarizedExperiment::rowData(reference)$ecdf,
        use_na=FALSE
      )[,1]
    )

  # Make the plot
  p<-ggplot2::ggplot(
    gdf,
    ggplot2::aes(x=assay_label, y=target, color=100*percentile, label=assay_label)) +

    ggplot2::geom_segment(
      ggplot2::aes(y = min, x = assay_label, yend = max,  xend = assay_label),
      color = "black"
    ) +
    # NB: This will take more effort, since the min/max is enough for the line,
    # but I'll need to compute from the whole population which isn't stored in
    # gdf currently.
    #ggridges::geom_density_ridges(aes(x=`Protein_Peptide`, y = (max-min))) +
    ggplot2::scale_color_gradient2(
      low="blue", high="red", midpoint = stats::median(100*gdf$percentile, na.rm=T)
    ) +
    #scale_color_gradient2() +

    ggplot2::labs(title=patient_name) +
    ggplot2::ylab(quote("Log2 amol/"*mu*"g")) +
    ggplot2::xlab("") +
    ggplot2::geom_point(
      ggplot2::aes(y=target),stat='identity', size=8, stroke=2, shape=21, fill="white"
    )  +
    ggplot2::geom_point(ggplot2::aes(y=target), stat="identity", size=8, shape=21, color="black") +
    ggplot2::geom_point(ggplot2::aes(y=target), stat="identity", size=10, shape=21, color="black") +
    # NB: This should be either the quantile or "ND" if 0.
    ggplot2::geom_text(
      ggplot2::aes(y=target, fontface="bold"),
      size=4,
      label=ifelse(gdf$target==0, "ND",format(round(100*gdf$percentile))),
      col="black"
    ) +
    # xlim(-10, 82) +
    ggplot2::scale_x_discrete(expand=ggplot2::expansion(mult=0.05)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none", axis.text.y = ggplot2::element_text(size=8, face="bold"))

  p
}

globalVariables(c("assay_label","target","percentile","min","max"))

