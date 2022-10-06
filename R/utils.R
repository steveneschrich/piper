clean_assay_name <- function(s) {
  stringr::str_replace_all(s, "_","\n")
}
