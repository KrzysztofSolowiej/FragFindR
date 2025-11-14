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

      mz_clicked <- as.numeric(click_data$x)
      if (is.na(mz_clicked)) return()

      mzs <- selected_mz()
      if (length(mzs) >= 2) mzs <- numeric()
      selected_mz(c(mzs, mz_clicked))
    })

    peaks <- eventReactive(input$custom_button, {
      req(input$custom_mz)

      txt <- input$custom_mz
      txt <- trimws(txt)
      if (nchar(txt) == 0) {
        validate(need(FALSE, "No data provided"))
      }

      mat <- tryCatch(
        utils::read.table(text = txt, header = FALSE, stringsAsFactors = FALSE, comment.char = "", blank.lines.skip = TRUE),
        error = function(e) NULL
      )

      validate(
        need(!is.null(mat) && ncol(mat) >= 2, "Input must contain at least two columns per line (m/z and intensity)")
      )

      mz <- as.numeric(mat[[1]])
      intensity <- as.numeric(mat[[2]])

      validate(
        need(!any(is.na(mz)), "m/z column contains non-numeric values"),
        need(!any(is.na(intensity)), "intensity column contains non-numeric values"),
        need(length(mz) >= 1, "No valid rows parsed")
      )

      data.frame(mz = mz, intensity = intensity)
    })

    output$custom_plot <- plotly::renderPlotly({
      req(peaks())
      df <- peaks()

      x_range <- range(df$mz, na.rm = TRUE)
      buf <- ifelse(diff(x_range) == 0, x_range[1] * 0.01 + 1e-6, 0.05 * diff(x_range))
      x_range <- c(x_range[1] - buf, x_range[2] + buf)

      p <- plotly::plot_ly(
        data = df,
        source = ns("custom_plot")
      ) %>%
        plotly::add_segments(
          x = ~mz, xend = ~mz,
          y = 0, yend = ~intensity,
          line = list(color = "#8F8F8F"),
          hoverinfo = "text",
          text = ~paste0("m/z: ", round(mz, 4), "<br>Intensity: ", round(intensity, 2))
        )

      if (identical(input$fix_labels, "Always")) {
        annotations <- lapply(seq_len(nrow(df)), function(i) {
          list(
            x = df$mz[i],
            y = df$intensity[i],
            text = formatC(df$mz[i], digits = 5, format = "f"),
            showarrow = FALSE,
            textangle = -90,
            xanchor = "center",
            yanchor = "bottom",
            font = list(size = 10)
          )
        })
        p <- p %>% plotly::layout(annotations = annotations)
      }

      p <- p %>% plotly::layout(
        title = NULL,
        xaxis = list(title = "m/z", range = x_range),
        yaxis = list(title = "Intensity"),
        showlegend = FALSE
      )

      plotly::event_register(p, "plotly_click")
    })

    output$custom_distance_info <- renderUI({
      mzs <- selected_mz()
      fmt <- function(x) ifelse(is.null(x), "—", sprintf("%.5f", x))
      first  <- if (length(mzs) >= 1) mzs[1] else NULL
      second <- if (length(mzs) >= 2) mzs[2] else NULL
      diff_v <- if (length(mzs) == 2) abs(diff(mzs)) else NULL

      div(
        style = "
          border: 1px solid #ddd;
          border-radius: 8px;
          padding: 8px 12px;
          background-color: #fafafa;
          font-size: 13px;
          line-height: 1.4;
          min-height: 80px;
        ",
        if (length(mzs) == 0) {
          tags$span("Click peaks to measure Δm/z between two fragments.",
                    style = 'color:#777; font-style:italic;')
        } else {
          tagList(
            div(tags$b("First m/z: "),  fmt(first)),
            div(tags$b("Second m/z: "), fmt(second)),
            div(tags$b("Δm/z: "),       ifelse(is.null(diff_v), "—", signif(diff_v, 6)))
          )
        }
      )
    })
  })
}
