#' Title
#'
#' @param .x
#' @param pnum
#'
#' @return
#' @export
#'
#' @importFrom stats median
#' @examples
plot_subway <- function(sample, reference, panels) {


  #panel_scores <- purrr::map(panels(ref), ~generate_signature(Biobase::exprs(ref), gene_set=which_targets(ref,panel=.x)))

  df <- tibble::tribble(
    ~x, ~y, ~panel, ~name, ~percentile,
    0,2, "Overview", "Biopsies",0,
    1,1, "Immunologist's Panel", "Immune Hot/Cold",0,
    2,1, "Immunologist's Panel", "Checkpoint Blockade", percentile(sample, reference, panels$`Immune Checkpoint`) |> median(),
    3,1, "Immunologist's Panel", "ADC Targets",0,
    4,1, "Immunologist's Panel", "Cancer Antigens", percentile(sample, reference, panels$`Cancer Antigens`) |> max(),
    1,2, "Oncologist's Panel", "RTK", percentile(sample, reference, panels$RTK) |> max(),
    2,2, "Oncologist's Panel", "KRAS G12C",0,
    3,2, "Oncologist's Panel", "MAPK", percentile(sample, reference, panels$MAPK) |> median(),
    4,2, "Oncologist's Panel", "CDK", percentile(sample, reference, panels$CDK) |> median(),
    5,2, "Oncologist's Panel", "PI3K/AKT/MTOR",percentile(sample, reference, panels$`PI3K/AKT/MTOR`) |> median(),
    1,3, "Molecular Pathologist's Panel", "QC", percentile(sample,reference, panels$`Tissue QC`) |> median(),
    2,3, "Molecular Pathologist's Panel", "Histology",0,
    3,3, "Molecular Pathologist's Panel", "EMT",percentile(sample,reference,panels$EMT) |> median(),
    4,3, "Molecular Pathologist's Panel", "Proliferation", percentile(sample, reference, panels$Proliferation) |> median(),
    5,3, "Molecular Pathologist's Panel", "Metabolism", percentile(sample, reference, panels$Metabolism) |> median()

    # Molecular Phenotyping (Path)
    # Cancer Signaling (Oncologist)
    # Tumor-Immune Microenvironment (Immunologist)
  ) |>
    dplyr::mutate(name = stringr::str_wrap(name, 10))

  df_lines <- tibble::tribble(
    ~x, ~y, ~xend, ~yend,~panel,
    0,2,1,1, "Immunologist's Panel",
    0,2,1,2, "Oncologist's Panel",
    0,2,1,3, "Molecular Pathologist's Panel",
    1,1,2,1,"Immunologist's Panel",
    2,1,3,1,"Immunologist's Panel",
    3,1,4,1,"Immunologist's Panel",
    1,2,2,2, "Oncologist's Panel",
    2,2,3,2, "Oncologist's Panel",
    3,2,4,2, "Oncologist's Panel",
    4,2,5,2, "Oncologist's Panel",
    1,3,2,3, "Molecular Pathologist's Panel",
    2,3,3,3, "Molecular Pathologist's Panel",
    3,3,4,3, "Molecular Pathologist's Panel",
    4,3,5,3, "Molecular Pathologist's Panel"
  )

  cpalette <- c(
    "Immunologist's Panel"="#bf9000",
    "Oncologist's Panel"="#548235",
    "Molecular Pathologist's Panel"="#2e75b6"
  )
  p<-ggplot2::ggplot(df, ggplot2::aes(x=x,y=y,col=panel, label = name, fill = percentile)) +
    ggplot2::geom_segment(ggplot2::aes(x=x,y=y,xend=xend,yend=yend, col=panel), data=df_lines, inherit.aes=FALSE,
                          size=4, alpha=0.3) +
    ggplot2::scale_color_manual(values = cpalette) +
    ggplot2::scale_fill_gradient(low="green",high="red") +
    ggplot2::geom_point(size=15, shape=21, stroke = 4) +
    #ggrepel::geom_text_repel(point.padding=2) +
    ggplot2::geom_text(nudge_y=0.4) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::expand_limits(x=c(-1,6), y=c(0,4))

  #pp<-plotly::ggplotly(p)
  p
}
globalVariables(c("x","y","panel","name","xend","yend"))
