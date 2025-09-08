#' neuloss UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_neuloss_ui <- function(id) {
  ns <- NS(id)
  tagList(
    DT::dataTableOutput(ns("results_table")),
    plotly::plotlyOutput(ns("spectrum_plot"), height = "480px")
  )
}

#' neuloss Sidebar Function
#'
#' @description Sidebar UI for a shiny neuloss Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_neuloss_sidebar <- function(id) {
  ns <- NS(id)

  tagList(
    #numericInput(ns("target_diff"), "Target Difference", value = 123.05, step = 0.01),
    textInput(
      ns("target_diff"),
      "Target Difference",
      value = "123.05"
    ),
    radioButtons(
      ns("polarity"),
      "Polarity",
      choices = c("Positive", "Negative"),
      selected = "Positive",
      inline = TRUE
    ),
    radioButtons(
      ns("tol_type"),
      "Tolerance Type",
      choices = c("Dalton", "PPM"),
      selected = "Dalton",
      inline = TRUE
    ),
    uiOutput(ns("tolerance_ui")),
    numericInput(ns("min_int"), "Min Relative Intensity (%)", value = 5, step = 1),
    selectInput(
      ns("collision_energy_level"),
      "Collision Energy Level",
      choices = c("Any" = "ALL", "low", "med", "high"),
      selected = "ALL"
    ),
    checkboxInput(ns("filter_voltage"), "Filter by Voltage Range", value = FALSE),
    uiOutput(ns("voltage_range_ui")),
    radioButtons(
      ns("plot_labels"),
      "Show labels",
      choices = c("On hover", "Always"),
      selected = "On hover",
      inline = TRUE
    ),
    actionButton(ns("run_search"), "Run search", class = "btn-custom"),
    uiOutput(ns("selected_structure"))
  )
}
