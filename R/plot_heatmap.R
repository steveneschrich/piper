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

    subtype = .x$subtype.y,
    laterality = .x$laterality_combined,
    # grade = lscc$Grade.or.Differentiation.Conformed...CAP,
    P53 = .x$TP53,
    annotation_name_gp = grid::gpar(fontsize= 8),
    col = list(
      "subtype" = c("Inflamed" = "red3", "Mixed"="purple", "Redox"="orange"),
      "laterality" = c("Left"="darkgrey","Right"="snow2"),
      # "grade"= c("1"="green","2"="lightgreen","3"="red3","X"="gray"),
      "P53"=c("mut"="red","trunc"="black","wt"="lightgreen")
    )

  )

  pt_split <- rep("reference",ncol(.x))
  pt_split[pnum]<-"target"

  #numna<-apply(Biobase::exprs(.x),1,function(x){length(which(is.na(x)))})
  #.x <- .x[numna < 70,]
  x<-t(scale(t(Biobase::exprs(.x))))
  h <- ComplexHeatmap::Heatmap(
    x,
    clustering_distance_rows = protein_distance,
    clustering_distance_columns = protein_distance,
    cluster_rows=FALSE,
    row_split = stringr::str_wrap(Biobase::fData(.x)$Subcategory,10),
    #column_split = lscc$subtype.y,
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
    #,
    #left_annotation = ComplexHeatmap::rowAnnotation(
    #  Category = ComplexHeatmap::anno_block(
    #    gp = grid::gpar(fill="lightgrey"),
    #    labels = stringr::str_wrap(Biobase::fData(lscc)$Subcategory,10),
    #    labels_gp = grid::gpar(col="white",fontsize=6)
    #  )

  )
  #decorate_heatmap_body("cases", {
  #  i = which(colnames(mat) == "1961")
  ##  x = i/ncol(mat)
  #  grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2, lty = 2))
  #  grid.text("Vaccine introduced", x, unit(1, "npc") + unit(5, "mm"))
  #})

  h
}
