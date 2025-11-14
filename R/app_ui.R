#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    fluidPage(
      #titlePanel("FragFindR HMDB Fragment Search"),
      titlePanel(title = span("FragFindR", class = "title-font")),
      sidebarLayout(
        sidebarPanel(uiOutput("dynamic_sidebar"), width = 3),
        mainPanel(
          tabsetPanel(id = "main_tabs",
            tabPanel("Fragment search", mod_frag_ui("frag")),
            tabPanel("Neutral loss search", mod_neuloss_ui("neuloss")),
            tabPanel("Custom peaks visualization", mod_custom_ui("custom"))
        )
        )
      )
    )
  )
}


#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd

golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "FragFindR"
    ),
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),

    # Spinner container (initially hidden)
    tags$div(
      id = "spinner-container",
      tags$img(id = "spinner-logo", src = "www/logo.png")
    ),

    # JavaScript to show/hide spinner
    tags$script(HTML("
      Shiny.addCustomMessageHandler('show_spinner', function(show) {
        var spinner = document.getElementById('spinner-container');
        spinner.style.display = show ? 'flex' : 'none';
      });
    "))
  )
}
