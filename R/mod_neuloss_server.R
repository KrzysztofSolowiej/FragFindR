#' neuloss Server Functions
#'
#' @noRd
mod_neuloss_server <- function(id, con, hmdb_name_map, hmdb_mass_map) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    results_data <- reactiveVal(NULL)

    # --- cache environment for spectra ---
    spectra_cache <- new.env(parent = emptyenv())

    output$tolerance_ui <- renderUI({
      if (input$tol_type == "Dalton") {
        numericInput(ns("tolerance"), "Tolerance (Da)", value = 0.01, step = 0.001)
      } else {
        numericInput(ns("tolerance"), "Tolerance (PPM)", value = 10, step = 1)
      }
    })

    output$voltage_range_ui <- renderUI({
      if (input$filter_voltage) {
        sliderInput(ns("voltage_range"), "Collision Energy Voltage Range",
                    min = 0, max = 120, value = c(0, 120), step = 1)
      }
    })

    observeEvent(input$run_search, {
      session$sendCustomMessage("show_spinner", TRUE)
      output$spectrum_plot <- plotly::renderPlotly(NULL)
      output$selected_structure <- renderUI(NULL)

      res <- local({
        tol_value <- input$tolerance
        if (input$tol_type == "PPM") tol_value <- (tol_value / 1e6) * input$target_diff

        tmp <- find_peak_diff(
          con,
          target_diff = input$target_diff,
          tolerance = tol_value,
          polarity = input$polarity,
          collision_energy_level = if (input$collision_energy_level != "ALL") input$collision_energy_level else NULL,
          voltage_range = if (isTRUE(input$filter_voltage)) input$voltage_range else NULL,
          min_rel_intensity = input$min_int,
          mass_map = hmdb_mass_map
        )

        if (nrow(tmp) == 0) {
          results_data(data.frame())
          output$results_table <- DT::renderDT(data.frame())
          return(NULL)
        }

        tmp$compound_name <- get_compound_name(tmp$hmdb_id, hmdb_name_map)
        tmp$structure     <- get_compound_smiles(tmp$hmdb_id, hmdb_name_map)
        tmp
      })
      gc()

      if (is.null(res) || nrow(res) == 0) {
        session$sendCustomMessage("show_spinner", FALSE)
        return()
      }

      res <- res %>%
        dplyr::arrange(precursor_type, diff_from_precursor) %>%
        dplyr::select(
          compound_name,
          precursor_type,
          diff_from_precursor,
          precursor_mz,
          frag_mz,
          analyzer_type,
          ionization_type,
          collision_energy_level,
          collision_energy_voltage,
          frag_intensity,
          frag_rel_intensity,
          hmdb_id,
          structure,
          file
        ) %>%
        dplyr::rename(
          `Compound Name` = compound_name,
          `Precursor Type` = precursor_type,
          `Difference from Precursor` = diff_from_precursor,
          `Precursor m/z` = precursor_mz,
          `Fragment m/z` = frag_mz,
          `Analyzer Type` = analyzer_type,
          `Ionization Type` = ionization_type,
          `Collision Energy Level` = collision_energy_level,
          `Collision Energy Voltage` = collision_energy_voltage,
          `Fragment Intensity` = frag_intensity,
          `Fragment Rel Intensity (%)` = frag_rel_intensity,
          `HMDB Id` = hmdb_id,
          `SMILES` = structure,
          `File Name` = file
        ) %>%
        dplyr::mutate(Index = dplyr::row_number())

      res$`HMDB Id Link` <- paste0(
        "<a href='https://www.hmdb.ca/metabolites/",
        res$`HMDB Id`, "' target='_blank'>", res$`HMDB Id`, "</a>"
      )

      results_data(res)

      res <- res %>%
        dplyr::select(
          Index,
          `Compound Name`,
          `HMDB Id Link`,
          `Precursor Type`,
          `Difference from Precursor`,
          `Precursor m/z`,
          `Fragment m/z`,
          `Analyzer Type`,
          `Ionization Type`,
          `Collision Energy Level`,
          `Collision Energy Voltage`,
          `Fragment Intensity`,
          `Fragment Rel Intensity (%)`
        )

      output$results_table <- DT::renderDT(
        res,
        escape = FALSE,
        extensions = c("ColReorder"),
        options = list(scrollX = TRUE, scrollY = "260px", paging = FALSE, dom = "Bfrtip", colReorder = TRUE),
        selection = "single",
        rownames = FALSE
      )

      session$sendCustomMessage("show_spinner", FALSE)
      gc()
    })

    plot_triggers <- reactive({ list(rows = input$results_table_rows_selected, labels = input$plot_labels) })

    observeEvent(plot_triggers(), {
      idx <- input$results_table_rows_selected
      if (length(idx) == 0) return()

      res <- results_data()
      row <- res[idx, ]
      precursor_mz <- row$`Precursor m/z`
      frag_mz <- row$`Fragment m/z`
      file_id <- row$`File Name`
      hmdb_id <- stringr::str_extract(file_id, "HMDB\\d+")
      compound_name <- get_compound_name(hmdb_id, hmdb_name_map)

      # --- Use cache ---
      if (exists(file_id, envir = spectra_cache)) {
        peaks <- get(file_id, envir = spectra_cache)
      } else {
        peaks <- DBI::dbGetQuery(con, glue::glue("
          SELECT mz, intensity
          FROM spectra
          WHERE file = '{file_id}'
        "))
        assign(file_id, peaks, envir = spectra_cache)
      }

      x_range <- range(peaks$mz, na.rm = TRUE)
      buffer <- 0.02 * diff(x_range)
      x_range <- c(x_range[1] - buffer, x_range[2] + buffer)

      p <- plotly::plot_ly() %>%
        plotly::add_segments(data = peaks, x = ~mz, xend = ~mz, y = ~0, yend = ~intensity,
                             line = list(color = "#8F8F8F"),
                             hoverinfo = "text", text = ~paste0("m/z: ", round(mz,4), "<br>Intensity: ", round(intensity,2))) %>%
        plotly::add_segments(data = peaks %>% dplyr::filter(abs(mz - precursor_mz) < 0.001),
                             x = ~mz, xend = ~mz, y = ~0, yend = ~intensity,
                             line = list(color = "#FF5454", width = 3),
                             hoverinfo = "text",
                             text = ~paste0("m/z: ", round(mz,4), "<br>Intensity: ", round(intensity,2))) %>%
        plotly::add_segments(data = peaks %>% dplyr::filter(abs(mz - frag_mz) < 0.001),
                             x = ~mz, xend = ~mz, y = ~0, yend = ~intensity,
                             line = list(color = "#2761F5", width = 3),
                             hoverinfo = "text",
                             text = ~paste0("m/z: ", round(mz,4), "<br>Intensity: ", round(intensity,2))) %>%
        plotly::layout(title = NULL,
                       xaxis = list(title = "m/z", range = x_range),
                       yaxis = list(title = "Intensity"),
                       showlegend = FALSE)

      if (input$plot_labels == "Always") {
        annotations <- lapply(1:nrow(peaks), function(i) {
          list(x = peaks$mz[i], y = peaks$intensity[i], text = round(peaks$mz[i],4),
               showarrow = FALSE, textangle = -90, xanchor = "center", yanchor = "bottom",
               font = list(size = 10, color = "black"))
        })
        p <- p %>% plotly::layout(annotations = annotations)
      }

      output$spectrum_plot <- plotly::renderPlotly(p)

      smi <- row$`SMILES`
      img_tag <- ""
      if (!is.na(smi) && nzchar(smi)) {
        img_tag <- tryCatch(smiles_to_img_tag(smi, compound_name, width = 250, height = 250), error = function(e) "")
      }

      output$selected_structure <- renderUI({ tagList(hr(), HTML(img_tag)) })
    })
  })
}
