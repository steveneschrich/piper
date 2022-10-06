#' Title
#'
#' @param .x
#' @param panel
#'
#'
#' @return
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @examples
panel_target_list <- function(.x, panel="Tissue QC") {
  stopifnot("ExpressionSet" %in% class(.x))

  Biobase::fData(.x) %>%
    dplyr::filter(.data$Subcategory %in% panel)
}

#' Title
#'
#' @param .x
#' @param panel
#'
#' @return
#' @export
#'
#' @examples
which_targets <- function(.x, panel="Tissue QC") {
  rownames(panel_target_list(.x, panel))

}
