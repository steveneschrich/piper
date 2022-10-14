#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
protein_distance <- function(x,y) {
  stopifnot(length(x)==length(y))
  n <- length(x)

  # nb, in PDS is number of blanks. Number missing, here calculated as the
  # number of NA's after adding the two vectors.
  nb <- length(which(is.na(x+y)))

  # Base case, if only one element is in common you cannot compute the distance,
  # so we take the distance as the number of elements (max distance possible).
  # Or the limit of Partial Distance Strategy (N/N-NB).
  if (nb < 2 ) return(n)

  # Otherwise, compute distance dropping observations
  d <- 1 - stats::cor(x, y, use="pairwise.complete.obs")

  # Weight d by missingness per PDS approach.
  d * (n/(n-nb))
}

#' Title
#'
#' @param .x
#' @param pnum
#' @param panel
#' @return
#' @export
#' @importFrom rlang .data
#' @examples
plot_heatmap<-function(.x, pnum=1, panel="Tissue QC") {

  .x <- .x[which_targets(.x, panel),]

  if ( is.character(pnum))
    pnum <- which(Biobase::sampleNames(.x) %in% pnum)


  ca <- ComplexHeatmap::columnAnnotation(

    subtype = .x$subtype,
    laterality = .x$laterality,
    annotation_name_gp = grid::gpar(fontsize= 8),
    col = list(
      "subtype" = c("Inflamed" = "red3", "Mixed"="purple", "Redox"="orange"),
      "laterality" = c("Left"="darkgrey","Right"="snow2")
    )

  )

  pt_split <- rep("reference",ncol(.x))
  pt_split[pnum]<-"target"

  x<-t(scale(t(Biobase::exprs(.x))))
  h <- ComplexHeatmap::Heatmap(
    x,
    clustering_distance_rows = protein_distance,
    clustering_distance_columns = protein_distance,
    cluster_rows=FALSE,
    row_split = stringr::str_wrap(Biobase::fData(.x)$Subcategory,10),
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    row_gap = grid::unit(2, "mm"),
    row_title_gp = grid::gpar(fontsize=8),
    row_title_rot = 0,
    top_annotation = ca,
    show_column_dend = FALSE,
    show_heatmap_legend = FALSE,
    show_row_names = TRUE,
    row_names_side="right",
    column_split =pt_split

  )

  h
}
