#' frag UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_frag_ui <- function(id) {
  ns <- NS(id)
  tagList(
    DT::dataTableOutput(ns("matches_table")),
    plotly::plotlyOutput(ns("spectrum_plot"), height = "480px")
  )
}

mod_frag_sidebar <- function(id) {
  ns <- NS(id)

  tagList(
    textInput(
      ns("fragments_mz"),
      "Fragment M/Z (separated by comma)",
      value = "184.07"
    ),
    radioButtons(
      ns("frag_polarity"),
      "Polarity",
      choices = c("Positive", "Negative"),
      selected = "Positive",
      inline = TRUE
    ),
    radioButtons(
      ns("frag_tol_type"),
      "Tolerance Type",
      choices = c("Dalton", "PPM"),
      selected = "Dalton",
      inline = TRUE
    ),
    uiOutput(ns("frag_tolerance_ui")),
    numericInput(ns("frag_min_int"), "Min Relative Intensity (%)", value = 5, step = 1),
    selectInput(
      ns("frag_collision_energy_level"),
      "Collision Energy Level",
      choices = c("Any" = "ALL", "low", "med", "high"),
      selected = "ALL"
    ),
    checkboxInput(ns("frag_filter_voltage"), "Filter by Voltage Range", value = FALSE),
    uiOutput(ns("frag_voltage_range_ui")),
    radioButtons(
      ns("frag_plot_labels"),
      "Show labels",
      choices = c("On hover", "Always"),
      selected = "On hover",
      inline = TRUE
    ),
    actionButton(ns("frag_run_search"), "Run search", class = "btn-custom"),
    uiOutput(ns("mz_distance_info")),
    uiOutput(ns("selected_structure"))
  )
}
