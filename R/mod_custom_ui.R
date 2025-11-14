#' custom UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_custom_ui <- function(id) {
  ns <- NS(id)
  tagList(
    plotly::plotlyOutput(ns("custom_plot"), height = "480px")
  )
}


mod_custom_sidebar <- function(id) {
  ns <- NS(id)

  tagList(
    textInput(
      ns("custom_mz"),
      "Custom m/z (separated by comma)",
      value = NULL
    ),
    textInput(
      ns("custom_intensity"),
      "Custom intensities (separated by comma)",
      value = NULL
    ),
    actionButton(ns("custom_button"), "Visualize", class = "btn-custom"),
    radioButtons(
      ns("fix_labels"),
      "Show labels",
      choices = c("On hover", "Always"),
      selected = "On hover",
      inline = TRUE
    ),
    uiOutput(ns("custom_distance_info"))
  )
}
