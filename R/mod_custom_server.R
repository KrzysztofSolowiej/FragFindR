#' custom Server Functions
#'
#' @noRd
mod_custom_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    selected_mz <- reactiveVal(numeric())

    observeEvent(plotly::event_data("plotly_click", source = ns("custom_plot")), {
      click_data <- plotly::event_data("plotly_click", source = ns("custom_plot"))
      if (is.null(click_data) || is.null(click_data$x)) return()

      mz_clicked <- click_data$x
      mzs <- selected_mz()
      if (length(mzs) >= 2) {
        mzs <- numeric()
      }
      mzs <- c(mzs, mz_clicked)
      selected_mz(mzs)
    })

    peaks <- eventReactive(input$custom_button, {
      req(input$custom_mz, input$custom_intensity)

      mz <- as.numeric(trimws(strsplit(input$custom_mz, ",")[[1]]))
      int <- as.numeric(trimws(strsplit(input$custom_intensity, ",")[[1]]))

      validate(
        need(length(mz) == length(int), "m/z and intensity must have the same length")
      )

      data.frame(mz = mz, intensity = int)
    })

    # Build plot reactively (updates on fix_labels change too)
    output$custom_plot <- plotly::renderPlotly({
      req(peaks())

      df <- peaks()

      p <- plotly::plot_ly(
        source = ns("custom_plot"),
        data = df
      ) %>%
        plotly::add_segments(
          x = ~mz, xend = ~mz,
          y = 0, yend = ~intensity,
          line = list(color = "#8F8F8F"),
          hoverinfo = "text",
          text = ~paste0(
            "m/z: ", round(mz, 4),
            "<br>Intensity: ", round(intensity, 2)
          )
        )

      # If labels always visible
      if (input$fix_labels == "Always") {
        annotations <- lapply(1:nrow(df), function(i) {
          list(
            x = df$mz[i],
            y = df$intensity[i],
            text = round(df$mz[i], 4),
            showarrow = FALSE,
            textangle = -90,
            xanchor = "center",
            yanchor = "bottom",
            font = list(size = 10, color = "black")
          )
        })

        p <- p %>% plotly::layout(annotations = annotations)
      }

      p %>% plotly::layout(
        xaxis = list(title = "m/z"),
        yaxis = list(title = "Intensity"),
        showlegend = FALSE
      )

    })

    output$custom_distance_info <- renderUI({
      mzs <- selected_mz()
      format_mz <- function(x) ifelse(is.null(x), "—", sprintf("%.5f", x))
      first <- if (length(mzs) >= 1) mzs[1] else NULL
      second <- if (length(mzs) >= 2) mzs[2] else NULL
      diff_val <- if (length(mzs) == 2) abs(diff(mzs)) else NULL

      tagList(
        br(),
        div(
          style = "
        border: 1px solid #ddd;
        border-radius: 8px;
        padding: 8px 12px;
        background-color: #fafafa;
        font-size: 13px;
        line-height: 1.4;
        min-height: 80px;
        display: flex;
        flex-direction: column;
        justify-content: center;
      ",
          if (length(mzs) == 0) {
            tags$span("Click peaks to measure Δm/z between two fragments.",
                      style = "color:#777; font-style:italic;")
          } else {
            tagList(
              tags$div(
                tags$b("First m/z: "), format_mz(first)
              ),
              tags$div(
                tags$b("Second m/z: "), format_mz(second)
              ),
              tags$div(
                tags$b("Δm/z: "), ifelse(is.null(diff_val), "—", signif(diff_val, 6))
              )
            )
          }
        )
      )
    })

  })
}
