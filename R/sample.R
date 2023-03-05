


#' Title
#'
#' @param .x
#' @param pnum
#' @param panel
#' @param use_na
#'
#' @return
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @examples
get_sample <- function(.x, pnum = 1, panel = NA, use_na = TRUE) {
  stopifnot(class(.x) %in% c("ExpressionSet"))

  # Extract specific sample (either the value or name, shouldn't matter)
  res <- Biobase::exprs(.x)[,pnum]%>%
    # Turn into tibble for filtering
    tibble::enframe(name = "Protein_Peptide", value = "Patient") %>%
    # Add assay annotation
    dplyr::left_join(
      Biobase::fData(.x),
      by = "Protein_Peptide"
    )

  # Filter by panel if requested
  if (!is.na(panel)) {
    res<- res %>% dplyr::filter(.data$Subcategory %in% panel)
  }

  # Some will be NA, use 0 in this instance (if requested)
  if ( ! use_na ) {
    res <- res %>% tidyr::replace_na(list(Patient = 0))

  }

  res
}

# How about make_sampleref_set(sample, ref). How much does it take to
# create it, calculate the statistics and go? Assuming ref is precomputed. You'd
# somehow need some extra stuff.

import_sample <- function(.x) {
  # Is it a vector of measurements with names?
  # Check nrow(.x) = nrow(reference)
  # Check rownames(.x) = rownames(reference)
  #
}

import_reference <- function(.x) {
  # No idea.
  # Calculate statistics, min, max

}
