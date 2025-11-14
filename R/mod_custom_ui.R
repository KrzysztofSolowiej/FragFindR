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
    textAreaInput(
      ns("custom_mz"),
      "Please paste your custom data (two columns: m/z intensity)",
      value = "90.97445 681
106.94476 274
110.02750 110
115.98965 95
117.98540 384
124.93547 613
124.99015 146
125.99793 207
133.95592 777
143.98846 478
144.99625 352
146.00410 999
151.94641 962
160.96668 387
163.00682 782
172.99055 17
178.95724 678
178.97725 391
180.97293 999
196.96778 720
208.96780 999
236.96245 999
254.97312 999",
      rows = 12, placeholder = "Paste two columns: m/z [whitespace] intensity, one pair per line"
    ),
    actionButton(ns("custom_button"), "Visualize", class = "btn-custom"),
    radioButtons(
      ns("fix_labels"),
      "Show labels",
      choices = c("On hover", "Always"),
      selected = "On hover",
      inline = TRUE
    ),
    br(),
    uiOutput(ns("custom_distance_info"))
  )
}
