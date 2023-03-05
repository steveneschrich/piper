
#' Determine if marker is in a list of panels
#'
#' @description A list of panels each contains a set of markers. Test
#' if a given marker is in each of the panels. Returns a (named) logical vector
#' indicating if the marker is the panel.
#'
#' @param marker A string representing a marker
#' @param panels A list of panels, each of which has a markers component.
#'
#' @return A list of named logicals, with the name representing the
#'  panel and the value representing if the marker is within the panel.
#' @export
#'
#' @examples
#' \dontrun{
#' markers_in_panels("CDN2A_ALLEAGALPNAPNSYGR", panels)
#' }
is_marker_in_panels <- function(marker, panels) {
  purrr::map_lgl(panels, \(p) {
    marker %in% p[["markers"]]
  })
}

#' Title
#'
#' @param marker
#' @param panels
#'
#' @return
#' @export
#'
#' @examples
marker_in_panels <- function(marker, panels) {
  in_panels <- is_marker_in_panels(marker, panels)
  names(panels)[in_panels] |>
    paste(collapse = ", ")
}

#' Title
#'
#' @param markers
#' @param panels
#'
#' @return
#' @export
#'
#' @examples
markers_in_panels <- function(markers, panels) {
  assertthat::assert_that(is.character(markers))
  assertthat::assert_that(is.list(panels))

  purrr::map_chr(markers, \(m) marker_in_panels(m, panels))
}

