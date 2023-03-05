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
  if (n-nb < 2 ) return(n)

  # Otherwise, compute distance dropping observations
  d <- 1 - stats::cor(x, y, use="pairwise.complete.obs")

  # Weight d by missingness per PDS approach.
  d * (n/(n-nb))
}

#' Title
#'
#' @param sample
#' @param reference
#' @param panel
#' @return
#' @export
#' @importFrom rlang .data
#' @examples
plot_heatmap<-function(sample, reference, panel) {
  #.x, pnum=1, panel="Tissue QC") {
  # Checking input parameters
  assertthat::assert_that(methods::is(reference, "SummarizedExperiment"))
  assertthat::assert_that(methods::is(sample, "SummarizedExperiment"))

  assertthat::assert_that(is.list(panel))
  assertthat::assert_that(utils::hasName(panel,"markers"))
  assertthat::assert_that(length(panel$markers)>0)

  # Panel: Select the specific panel to plot
  markers <- match_markers(panel$markers, reference)

  # Subset the sample/reference data
  sample <- sample[markers,]
  reference <- reference[markers,]

  ca <- ComplexHeatmap::columnAnnotation(
    subtype = c(sample$subtype, reference$subtype),
    laterality = c(sample$laterality, reference$laterality),
    annotation_name_gp = grid::gpar(fontsize= 8),
    col = list(
      "subtype" = c("Inflamed" = "red3", "Mixed"="purple", "Redox"="orange"),
      "laterality" = c("Left"="darkgrey","Right"="snow2")
    )

  )

  pt_split <- c("target",rep("reference",ncol(reference)))

  x <- cbind(
    SummarizedExperiment::assay(sample),
    SummarizedExperiment::assay(reference)
  )

  x<-t(scale(t(x)))
  h <- ComplexHeatmap::Heatmap(
    x,
    clustering_distance_rows = protein_distance,
    clustering_distance_columns = protein_distance,
    cluster_rows=FALSE,
    row_split = stringr::str_wrap(rep(panel$name,nrow(x)),10),
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
