clean_assay_name <- function(s) {
  stringr::str_replace_all(s, "_","\n")
}


#' Title
#'
#' @param marker_list
#' @param reference
#'
#' @return
#' @export
#'
#' @examples
match_markers <- function(markers, reference) {

  assertthat::assert_that(is.vector(markers))

  # The markers can be one of two things: either a numeric list of indices or
  # a character list of names. If numeric, we have to trust the caller knew what
  # they were doing. If character, we find the intersection.
  if ( all(is.numeric(markers)) ) {
    assertthat::assert_that(length(markers) <= nrow(reference))
    markers <- rownames(reference)[markers]
  }

  markers <- intersect(markers, rownames(reference))
  assertthat::assert_that(length(markers)>0, msg = "No panel markers identified in reference data set.")


  markers
}
