
#' Plot assay results using radar plot
#'
#' @description Plot an assay panel result for a sample against a reference
#'  set using a radar plot.
#'
#' @param sample A SummarizedExperiment object for the sample
#' @param reference A SummarizedExperiment object representing the reference population.
#' @param panel A panel list (markers, name, category)
#'
#' @return A [plotly::plot_ly()] radar plot representing the specified
#' panel for a target sample against a reference population.
#'
#' @export
#'
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' plot_radar(lscc[,1],lscc,
#'   panel = list(name="ADC",markers=list("PSCA_AVG","TACD2_AAG"),
#'   category = "Cell Surface, RTKs & ADC Targets"))
#' }
plot_radar <- function(sample, reference, panel=NULL,
                       axis.label.size=4) {

  # Checking input parameters
  assertthat::assert_that(methods::is(reference, "SummarizedExperiment"))
  assertthat::assert_that(methods::is(sample, "SummarizedExperiment"))
  assertthat::assert_that(nrow(sample) == nrow(reference))

  # Handle the panel input
  if (is.null(panel))
    panel <- list(
      "markers"=sample(rownames(reference),5),
      category = "Random", name = "Random"
      )
  assertthat::assert_that(is.list(panel))
  assertthat::assert_that(utils::hasName(panel,"markers"))
  assertthat::assert_that(length(panel$markers)>0)
  assertthat::assert_that(all(panel$markers %in% rownames(reference)))

  # Panel: Select the specific panel to plot
  markers <- match_markers(panel$markers, reference)
  panel_name <- panel$name

  # Subset the sample/reference data
  sample <- sample[markers,]
  reference <- reference[markers,]

  # Statistics: We plot the percent of max abundance in the radar plot.
  SummarizedExperiment::assays(sample)$percent_max <- row_percent(
    sample,
    ecdfs = SummarizedExperiment::rowData(reference)$ecdf,
    use_na=FALSE
  )

  qc <- tibble::tibble(
    assay_id = rownames(sample),
    split_assay_id(.data$assay_id),
    assay_label = label_format(.data$assay_id),
    percentile = row_percent(
      sample,
      ecdfs = SummarizedExperiment::rowData(reference)$ecdf,
      use_na=FALSE
    )[,1],
    abundance = SummarizedExperiment::assays(sample)$exprs[,1],
    text = sprintf("Protein: %s\nPeptide: %s\nAbundance: %5.2f\nPercentile: %5.2f %%",
                   .data$Protein, .data$Peptide,
                   .data$abundance, 100*.data$percentile)
  )

  # NB: Add an extra row with the first element in order to get points to connect. The
  # repeated row seems to allow that.
  qc <- dplyr::bind_rows(qc, dplyr::slice(qc, 1))

  p <- radar_init_plot(qc, title = paste(panel$name,panel$category,sep="\n")) |>
    radar_add_assay(q=~percentile, assay=~assay_label, text = ~text)

  p
}


#' Title
#'
#' @param p
#'
#' @return
#'
#' @examples
radar_add_assay <- function(p, q, assay, text) {
  plotly::add_trace(
    p,
    r = q, # Precalculated percentage
    theta = assay, # Precalculated angles
    text = text, # hover labels
    hoverinfo = "text",
    mode = "lines+markers", # We want both points and lines connected points, leave text for hover
    showlegend = FALSE, # No legend
    #cliponaxis=FALSE,
    connectgaps = TRUE,

    # Marker (point) characteristics
    marker = list(
      size = 20,
      color = "red"
    ),
    # Connector line characteristics
    line = list(
      color = "red",
      width = 5
    )
  )
}

#' Title
#'
#' @param p
#'
#' @return
#'
#' @examples
radar_init_plot<-function(x, title = "") {
  p <- plotly::plot_ly(x, type = "scatterpolar", mode = "none")

  plotly::layout(
    p,
    title = list(
      text = title,
      x = 0
    ),
    autosize = TRUE,
    margin =list(t=80,pad=1,b=80),
    polar = list(
      bgcolor = "#F5F5F4", # Light gray
      # Radial axis are the concentric circles
      radialaxis = list(
        range = c(-0.1,1.01), # Slightly larger to accomodate tickmark at 0, 1.0 with space.
        tickvals = list(0,0.5,1), # Just ticks at center, halfway and outer
        ticktext = list("0%","50%","100%"), # Label in percentage rather than decimal
        griddash = I("dash"), # Doesn't seem to work
        ticks = "", # Don't draw tickmarks
        tickfont = list(size = 14) # Change default font size of labels (0,50,100%).
      ),
      # These are the spokes
      angularaxis = list(
        gridcolor = "lightgray", # Lightgray spokes coming from origin
        griddash = I("dash"), # Doesn't seem to work
        gridwidth = 0,
        tickfont = list(size =15),
        #showticklabels = TRUE, # Turn off tick labels (just angles, not needed).
        ticks = "", # We don't want the tick, in addition to tick labels
        showline = FALSE # Turns off outer ring (handle via a radialaxis tickmark)
      )
    )
  )

}
