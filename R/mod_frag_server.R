#' frag Server Functions
#'
#' @noRd
mod_frag_server <- function(id, con, hmdb_name_map, hmdb_mass_map) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Store results after search
    results_data <- reactiveVal(NULL)

    output$frag_tolerance_ui <- renderUI({
      if (input$frag_tol_type == "Dalton") {
        numericInput(ns("frag_tolerance"), "Tolerance (Da)", value = 0.01, step = 0.001)
      } else {
        numericInput(ns("frag_tolerance"), "Tolerance (PPM)", value = 10, step = 1)
      }
    })

    output$frag_voltage_range_ui <- renderUI({
      if (input$frag_filter_voltage) {
        sliderInput(
          ns("frag_voltage_range"),
          "Collision Energy Voltage Range",
          min = 0, max = 120, value = c(0, 120), step = 1
        )
      }
    })

    # Run search when button clicked
    observeEvent(input$frag_run_search, {
      output$spectrum_plot <- plotly::renderPlotly(NULL)
      output$selected_structure <- renderUI(NULL)

      # Parse target m/z values to compute ppm if needed
      frag_vals <- strsplit(input$fragments_mz, ",")[[1]] %>%
        trimws() %>% as.numeric()

      tol_value <- input$frag_tolerance

      if (input$frag_tol_type == "PPM") {
        # Use mean of fragment m/z values for conversion
        tol_value <- (tol_value / 1e6) * mean(frag_vals, na.rm = TRUE)
      }

      matches <- if (input$frag_tol_type == "PPM") {
        find_fragments(con,
                       fragments_mz = input$fragments_mz,
                       tolerance_ppm = input$frag_tolerance,
                       polarity = input$frag_polarity,
                       collision_energy_level = if (input$frag_collision_energy_level != "ALL") input$frag_collision_energy_level else NULL,
                       voltage_range = if (isTRUE(input$frag_filter_voltage)) input$frag_voltage_range else NULL,
                       min_rel_intensity = input$frag_min_int)
      } else {
        find_fragments(con,
                       fragments_mz = input$fragments_mz,
                       tolerance_da = input$frag_tolerance,
                       polarity = input$frag_polarity,
                       collision_energy_level = if (input$frag_collision_energy_level != "ALL") input$frag_collision_energy_level else NULL,
                       voltage_range = if (isTRUE(input$frag_filter_voltage)) input$frag_voltage_range else NULL,
                       min_rel_intensity = input$frag_min_int)
      }

      if (nrow(matches) == 0) {
        results_data(data.frame())  # store empty df
        output$matches_table <- DT::renderDT(data.frame())
        output$spectrum_plot <- plotly::renderPlotly(NULL)
        output$selected_structure <- renderUI(NULL)
        return()
      }

      # matches$compound_name <- get_compound_name(matches$hmdb_id, hmdb_name_map)
      # matches$structure <- get_compound_smiles(matches$hmdb_id, hmdb_name_map)
      # matches$structure_img <- vapply(
      #   seq_len(nrow(matches)),
      #   function(i) {
      #     smi <- matches$structure[i]
      #     nm  <- matches$compound_name[i]
      #     if (!is.na(smi) && nzchar(smi)) {
      #       tryCatch(
      #         smiles_to_img_tag(smi, nm, width = 150, height = 150),
      #         error = function(e) ""
      #       )
      #     } else {
      #       ""
      #     }
      #   },
      #   FUN.VALUE = character(1)
      # )

      matches$compound_name <- get_compound_name(matches$hmdb_id, hmdb_name_map)
      matches$structure <- get_compound_smiles(matches$hmdb_id, hmdb_name_map)

      matches$structure_img <- vapply(
        seq_len(nrow(matches)),
        function(i) {
          smi <- matches$structure[i]
          nm  <- matches$compound_name[i]
          if (!is.na(smi) && nzchar(smi)) {
            tryCatch(smiles_to_img_tag(smi, nm, width = 150, height = 150), error = function(e) "")
          } else ""
        },
        FUN.VALUE = character(1)
      )

      #hmdb_mass_map <- readRDS("inst/extdata/hmdb_mass_map.rds")
      matches <- matches %>%
        dplyr::left_join(hmdb_mass_map %>%
                  dplyr::select(database_id, adduct_mass),
                  by = c("hmdb_id" = "database_id"))

      matches <- matches %>%
        dplyr::select(
          compound_name,
          match_count,
          hmdb_id,
          polarity,
          adduct_mass,
          analyzer_type,
          ionization_type,
          collision_energy_level,
          collision_energy_voltage,
          structure,
          file
        )       %>%
        dplyr::rename(
          `Compound Name` = compound_name,
          `Match Count` = match_count,
          `HMDB Id` = hmdb_id,
          `Polarity` = polarity,
          `Precursor m/z` = adduct_mass,
          `Analyzer Type` = analyzer_type,
          `Ionization Type` = ionization_type,
          `Collision Energy Level` = collision_energy_level,
          `Collision Energy Voltage` = collision_energy_voltage,
          `SMILES` = structure,
          `File Name` = file
        )

      matches <- matches %>%
        dplyr::mutate(Index = as.numeric(rownames(matches))) %>%
        dplyr::select(Index, everything())

      matches$`HMDB Id Link` <- paste0(
        "<a href='https://www.hmdb.ca/metabolites/",
        matches$`HMDB Id`,
        "' target='_blank'>",
        matches$`HMDB Id`,
        "</a>"
      )

      matches <- matches %>%
        dplyr::arrange(`Compound Name`) %>%
        dplyr::mutate(Index = dplyr::row_number())

      results_data(matches)

      matches <- matches %>%
        dplyr::select(
          Index,
          `Compound Name`,
          `HMDB Id Link`,
          `Polarity`,
          `Precursor m/z`,
          `Analyzer Type`,
          `Ionization Type`,
          `Collision Energy Level`,
          `Collision Energy Voltage`
          )

      output$matches_table <- DT::renderDT(
        matches,
        escape = FALSE,
        extensions = c("ColReorder"),
        options = list(
          scrollX = TRUE,
          scrollY = "260px",
          paging = FALSE,
          dom = "Bfrtip",
          colReorder = TRUE
        ),
        selection = "single",
        rownames = FALSE
      )
    })

    frag_plot_triggers <- reactive({
      list(
        rows = input$matches_table_rows_selected,
        labels = input$frag_plot_labels
      )
    })

    observeEvent(frag_plot_triggers(), {
      idx <- input$matches_table_rows_selected
      if (length(idx) == 0) return()

      matches <- results_data()
      if (is.null(matches) || nrow(matches) == 0) {
        warning("No matches available.")
        return()
      }

      row <- matches[idx, , drop = FALSE]
      if (nrow(row) == 0 || is.na(row$`File Name`)) {
        warning("Selected row has no valid file_id.")
        return()
      }

      file_id <- row$`File Name`[1]  # enforce single row
      hmdb_id <- stringr::str_extract(file_id, "HMDB\\d+")
      compound_name <- get_compound_name(hmdb_id, hmdb_name_map)

      query <- glue::glue("SELECT mz, intensity FROM spectra WHERE file = '{file_id}'")

      peaks <- tryCatch({
        DBI::dbGetQuery(con, query)
      }, error = function(e) {
        warning("Query failed: ", conditionMessage(e))
        NULL
      })

      if (nrow(peaks) == 0) {
        warning("No peaks found for this file.")
        return()
      }

      x_range <- range(peaks$mz, na.rm = TRUE)
      buffer <- 0.02 * diff(x_range)
      x_range <- c(x_range[1] - buffer, x_range[2] + buffer)

      p <- plotly::plot_ly() %>%
        plotly::add_segments(
          data = peaks,
          x = ~mz, xend = ~mz,
          y = ~0, yend = ~intensity,
          line = list(color = "#8F8F8F"),
          hoverinfo = "text",
          text = ~paste0("m/z: ", round(mz, 4), "<br>Intensity: ", round(intensity, 2))
        )

      frag_vals <- strsplit(input$fragments_mz, ",")[[1]] |> trimws() |> as.numeric()

      tol_da_for <- function(mz0) {
        if (identical(input$frag_tol_type, "PPM")) {
          (input$frag_tolerance / 1e6) * mz0
        } else {
          input$frag_tolerance
        }
      }

      for (mz_val in frag_vals) {
        tol_this <- tol_da_for(mz_val)
        p <- p %>%
          plotly::add_segments(
            data = peaks %>% dplyr::filter(abs(mz - mz_val) <= tol_this),
            x = ~mz, xend = ~mz,
            y = ~0, yend = ~intensity,
            line = list(color = "#FF5454", width = 3),
            hoverinfo = "text",
            text = ~paste0(
              "m/z: ", round(mz, 4),
              "<br>Intensity: ", round(intensity, 2),
              "<br>|Δ|: ", signif(abs(mz - mz_val), 6), " Da",
              "<br>tol: ", signif(tol_this, 6), " Da"
            )
          )
      }

      if (input$frag_plot_labels == "Always") {
        annotations <- lapply(1:nrow(peaks), function(i) {
          list(
            x = peaks$mz[i],
            y = peaks$intensity[i],
            text = round(peaks$mz[i], 4),
            showarrow = FALSE,
            textangle = -90,
            xanchor = "center",
            yanchor = "bottom",
            font = list(size = 10, color = "black")
          )
        })

        p <- p %>%
          plotly::layout(annotations = annotations)
      }

      p <- p %>%
        plotly::layout(
          #title = paste0("MS/MS Spectrum: ", compound_name, " (", hmdb_id, ")"),
          title = NULL,
          xaxis = list(title = "m/z", range = x_range),  # fix axis range
          yaxis = list(title = "Intensity"),
          showlegend = FALSE
        )

      output$spectrum_plot <- plotly::renderPlotly(p)

      smi <- row$`SMILES`
      img_tag <- if (!is.na(smi) && nzchar(smi)) {
        tryCatch(smiles_to_img_tag(smi, compound_name, width = 250, height = 250), error = function(e) "")
      } else ""
      output$selected_structure <- renderUI({ tagList(hr(), HTML(img_tag)) })
    })

  })
}
