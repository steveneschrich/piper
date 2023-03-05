
#' Return markers present in the experiment
#'
#' @description Given a matrix of assay data, determine the intersection
#'  of markers and the assay data. Note that a numeric list can be
#'  used for markers, or the actual names.
#'
#' @param markers A list of markers, some of which may be in the
#'  experiment object.
#' @param reference A matrix of assay measurements with rownames that
#' may include the marker list.
#'
#' @return The markers that are present in the reference.
#' @export
#'
#' @examples
#' \dontrun{
#' match_markers(c("Sepal.Length","NotThere"), t(iris))
#' [1] "Sepal.Length"
#' }
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




#' Pad labels to max width
#'
#' @description Given input parameters (labels), find the largest label and pad
#'  all strings to that width. Then concatenate into a single string.
#'
#' @details When plotting figures in [plotly::plot_ly()] the padding
#' is not always great. This function assumes the goal is to create
#' a multi-line string consisting of elements of the input. It determines
#' the maximum width of input labels, then centers all labels to that width
#' using spaces. The argument `extra_width` can be used to increase the
#' padding beyond the widest string.
#'
#' The padded strings are then concatenated together into a single
#' string, with a newline separating the strings. Additionally, an
#' extra line is added before and after the strings.
#'
#' @param ... A series of labels (strings) to concatenate with the same width
#' @param extra_width Additional width to add to longest string component (default = 6)
#'
#' @return A string with extra padding and concatenation
#' @export
#'
#' @examples
#' \dontrun{
#' label_pad("Line1","This is Line2", extra_width=0)
#'
#' [1] "             \n    Line1    \nThis is Line2\n             "
#' }
label_pad <- function( ..., extra_width=12) {
  s <- list(...)
  max_str <- max(stringi::stri_length(s))
  stringi::stri_pad_both(c("",s,""), width = max_str + extra_width) |>
    stringi::stri_flatten(collapse="\n")
}


#' Format labels with padding on all sides
#'
#' @description The assay label is a composite field. This function splits the
#'  composite label into components and pads them to a max width. Then combines
#'  them into a single string.
#'
#' @details This function is a wrapper function for ease of use. It uses
#' [split_assay_id()] to split a composite label into
#'  components. Then it uses [label_pad()] to pad and concatenate the results into
#'  single strings.
#' @param s A vector of strings
#'
#' @return A vector of strings consisting of padded, line-separated components
#'  of the input assay fields.
#' @export
#'
#' @examples
#' \dontrun{
#' label_format("K2C5_LAELEEALQK")
#' [1] "                      \n         K2C5         \n      LAELEEALQK      \n                      "
#' }
label_format <- function(s) {
  split_assay_id(s) |>
    dplyr::mutate(label = purrr::pmap_chr(list(Protein, Peptide), label_pad)) |>
    dplyr::pull(label)

}

#' Split assay ID into protein and peptide
#'
#' @description The assay id is a protein and peptide sequence combined. This
#'  function splits them into a data frame of the two components.
#'
#' @param s A vector of assay id's.
#'
#' @return A data frame of `Protein` and `Peptide` columns corresponding
#' to the content of the input vector.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' split_assay_id("K2C7_SAYGGPVGAGIR")
#' }
split_assay_id <- function(s) {
  s <- stringi::stri_split_fixed(s, "_", 2, simplify = TRUE)
  colnames(s) <- c("Protein","Peptide")

  tibble::as_tibble(s)
}


