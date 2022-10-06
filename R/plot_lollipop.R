#' Title
#'
#' @param .x
#' @param pnum
#' @param panel
#'
#' @return
#' @export
#'
#' @examples
plot_lollipop <- function(.x, pnum = 1, panel = "Tissue QC") {
  if ( is.character(pnum))
    pnum <- which(Biobase::sampleNames(.x) %in% pnum)

  patient_name <- Biobase::sampleNames(.x)[pnum]
  x <- Biobase::exprs(.x) |>
    as.data.frame() |>
    tibble::rownames_to_column("Protein_Peptide") |>
    tidyr::pivot_longer(
      c(dplyr::everything(),-.data$Protein_Peptide),
      names_to="sample",
      values_to="expression"
    ) |>
    dplyr::group_by(.data$Protein_Peptide) |>
    dplyr::summarize(
      max = max(.data$expression, na.rm=T),
      min = min(.data$expression, na.rm=T),
      ecdf_fun = list(stats::ecdf(.data$expression))
    )


  rev_label <- function(s) {
    x <- stringr::str_split_fixed(s, "_",n= 2)
    paste0(x[,2],"_",x[,1])
  }

  sample <- Biobase::exprs(.x)[,pnum] |>
    tibble::enframe(name="Protein_Peptide", value="Patient") |>
    dplyr::left_join(Biobase::fData(.x),
                     by=c("Protein_Peptide"="Protein_Peptide"))

  # The plotting data is the combination of x and the
  # sample data, with the "sample_quantile" combining
  # data from both.
  gdf <- sample |>
    dplyr::left_join(x, by=c("Protein_Peptide"="Protein_Peptide")) |>
    dplyr::mutate(
      sample_quantile =
        purrr::map2_dbl(.data$ecdf_fun, .data$Patient, ~.x(.y)
        )
    ) |>
    dplyr::mutate(
      protein_label = rev_label(.data$Protein_Peptide),
      Patient = ifelse(is.na(.data$Patient), 0, .data$Patient)
    ) |>
    # Note sure why this is here
    #dplyr::filter(Subcategory %in% c("Cancer Antigens","CDK","MAPK","RTK","PI3K/AKT/MTOR","Tissue QC")) |>
    dplyr::mutate(Protein_Peptide=stringr::str_replace_all(.data$Protein_Peptide, "_","\n"))

  # Tmp filter

  gdf <- gdf |> dplyr::filter(.data$Subcategory==panel)
  #dplyr::slice(c(3,5,4,7,6,2,1)) %>%
  # dplyr::mutate(Protein_Peptide = forcats::fct_inorder(Protein_Peptide))

  # Make the plot
  p<-ggplot2::ggplot(gdf,
                     ggplot2::aes(x=`Protein_Peptide`, y=Patient, color=100*sample_quantile, label=Protein_Peptide)) +

    ggplot2::geom_segment(
      ggplot2::aes(y = min, x = Protein_Peptide, yend = max,  xend = Protein_Peptide),
      color = "black"
    ) +
    # NB: This will take more effort, since the min/max is enough for the line,
    # but I'll need to compute from the whole population which isn't stored in
    # gdf currently.
    #ggridges::geom_density_ridges(aes(x=`Protein_Peptide`, y = (max-min))) +
    ggplot2::scale_color_gradient2(
      low="blue", high="red", midpoint = stats::median(100*gdf$sample_quantile, na.rm=T)
    ) +
    #scale_color_gradient2() +

    ggplot2::labs(title=patient_name) +
    ggplot2::ylab(quote("Log2 amol/"*mu*"g")) +
    ggplot2::geom_point(
      ggplot2::aes(y=Patient),stat='identity', size=8, stroke=2, shape=21, fill="white"
    )  +
    ggplot2::geom_point(ggplot2::aes(y=Patient), stat="identity", size=8, shape=21, color="black") +
    ggplot2::geom_point(ggplot2::aes(y=Patient), stat="identity", size=10, shape=21, color="black") +
    # NB: This should be either the quantile or "ND" if 0.
    ggplot2::geom_text(
      ggplot2::aes(y=Patient, fontface="bold"),
      size=4,
      label=ifelse(gdf$Patient==0, "ND",format(round(100*gdf$sample_quantile))),
      col="black"
    ) +
    # xlim(-10, 82) +
    ggplot2::scale_x_discrete(expand=ggplot2::expansion(mult=0.05)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none", axis.text.y = ggplot2::element_text(size=8, face="bold"))

  p
}

globalVariables(c("Protein_Peptide","Patient","sample_quantile","min","max"))

