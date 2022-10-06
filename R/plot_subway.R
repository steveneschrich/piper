#' Title
#'
#' @param .x
#' @param pnum
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @examples
plot_subway <- function(.x, pnum=1) {

  ref <- .x
  #panel_scores <- purrr::map(panels(ref), ~generate_signature(Biobase::exprs(ref), gene_set=which_targets(ref,panel=.x)))

  df <- tibble::tribble(
    ~x, ~y, ~panel, ~name,
    0,2, "Overview", "Biopsies",
    1,1, "Immunologist's Panel", "Immune Hot/Cold",
    2,1, "Immunologist's Panel", "Checkpoint Blockade",
    3,1, "Immunologist's Panel", "ADC Targets",
    4,1, "Immunologist's Panel", "Cancer Antigens",
    1,2, "Oncologist's Panel", "RTK",
    2,2, "Oncologist's Panel", "KRAS G12C",
    3,2, "Oncologist's Panel", "MAPK",
    4,2, "Oncologist's Panel", "CDK",
    5,2, "Oncologist's Panel", "PI3K/AKT/MTOR",
    1,3, "Molecular Pathologist's Panel", "QC",
    2,3, "Molecular Pathologist's Panel", "Histology",
    3,3, "Molecular Pathologist's Panel", "EMT",
    4,3, "Molecular Pathologist's Panel", "Proliferation",
    5,3, "Molecular Pathologist's Panel", "Metabolism"

    # Molecular Phenotyping (Path)
    # Cancer Signaling (Oncologist)
    # Tumor-Immune Microenvironment (Immunologist)
  ) %>%
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
  p<-ggplot2::ggplot(df, ggplot2::aes(x=x,y=y,col=panel, label = name)) +
    ggplot2::geom_segment(ggplot2::aes(x=x,y=y,xend=xend,yend=yend, col=panel), data=df_lines, inherit.aes=FALSE,
                          size=4, alpha=0.3) +
    ggplot2::scale_color_manual(values = cpalette) +
    ggplot2::geom_point(size=15, shape=21, stroke = 4, fill="white") +
    #ggrepel::geom_text_repel(point.padding=2) +
    ggplot2::geom_text(nudge_y=0.4) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::expand_limits(x=c(-1,6), y=c(0,4))

  #pp<-plotly::ggplotly(p)
  p
}
globalVariables(c("x","y","panel","name","xend","yend"))
