#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
#' @importFrom magrittr %>%
app_server <- function(input, output, session) {

  hmdb_name_path <- system.file("extdata", "hmdb_names_smiles_map.rds", package = "FragFindR")
  hmdb_mass_path <- system.file("extdata", "hmdb_mass_map.rds", package = "FragFindR")

  if (hmdb_name_path == "" || !file.exists(hmdb_name_path)) {
    stop("hmdb_names_smiles_map.rds not found!")
  }
  if (hmdb_mass_path == "" || !file.exists(hmdb_mass_path)) {
    stop("hmdb_mass_map.rds not found!")
  }

  hmdb_name_map <- readRDS(hmdb_name_path)
  hmdb_mass_map <- readRDS(hmdb_mass_path)

  # Writable location for the DuckDB file
  db_path <- file.path(rappdirs::user_data_dir("FragFindR"), "hmdb_spectra.duckdb")

  if (!file.exists(db_path)) {
    message("Downloading database file...")
    dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
    download.file(
      url = "https://www.dropbox.com/scl/fi/y00om5xwdis5l5c3bi5fx/hmdb_spectra.duckdb?rlkey=oijk6jq9wve0txw564jorb91r&dl=1",
      destfile = db_path,
      mode = "wb"
    )
  }

  con <- DBI::dbConnect(duckdb::duckdb(), db_path, read_only = TRUE)

  # Call modules
  mod_frag_server("frag", con, hmdb_name_map, hmdb_mass_map)
  mod_neuloss_server("neuloss", con, hmdb_name_map, hmdb_mass_map)

  output$dynamic_sidebar <- renderUI({
    switch(input$main_tabs,
           "Fragment search" = mod_frag_sidebar("frag"),
           "Neutral loss search" = mod_neuloss_sidebar("neuloss")
    )
  })

  # Close connection on session end
  session$onSessionEnded(function() {
    DBI::dbDisconnect(con, shutdown = TRUE)
  })
}
