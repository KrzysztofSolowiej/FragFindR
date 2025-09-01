#' helpers
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd
find_fragments <- function(con,
                           fragments_mz,
                           tolerance_da = NULL,          # pass one of these
                           tolerance_ppm = NULL,         # (not both)
                           polarity = NULL,
                           analyzer_type = NULL,
                           ionization_type = NULL,
                           collision_energy_level = NULL,
                           voltage_range = NULL,
                           min_rel_intensity = 5) {

  # --- validate inputs ---
  if (is.null(tolerance_da) && is.null(tolerance_ppm)) {
    stop("Provide either tolerance_da or tolerance_ppm.")
  }
  if (!is.null(tolerance_da) && !is.null(tolerance_ppm)) {
    stop("Use only one of tolerance_da or tolerance_ppm (not both).")
  }

  target_mz <- fragments_mz |>
    strsplit(",") |> unlist() |> trimws() |> as.numeric()
  if (!length(target_mz) || anyNA(target_mz)) {
    stop("Invalid fragment m/z values.")
  }

  # --- build dynamic filters ---
  esc <- function(x) gsub("'", "''", x)  # simple SQL escape for single quotes
  filters <- c()
  if (!is.null(polarity) && nzchar(polarity))
    filters <- c(filters, glue::glue("LOWER(polarity) = LOWER('{esc(polarity)}')"))
  if (!is.null(analyzer_type) && nzchar(analyzer_type))
    filters <- c(filters, glue::glue("LOWER(analyzer_type) = LOWER('{esc(analyzer_type)}')"))
  if (!is.null(ionization_type) && nzchar(ionization_type))
    filters <- c(filters, glue::glue("LOWER(ionization_type) = LOWER('{esc(ionization_type)}')"))
  if (!is.null(collision_energy_level) && nzchar(collision_energy_level))
    filters <- c(filters, glue::glue("collision_energy_level = '{esc(collision_energy_level)}'"))
  if (!is.null(voltage_range) && length(voltage_range) == 2)
    filters <- c(filters, glue::glue("collision_energy_voltage BETWEEN {voltage_range[1]} AND {voltage_range[2]}"))

  where_clause <- if (length(filters)) paste("WHERE", paste(filters, collapse = " AND ")) else ""

  # --- tolerance condition (per-target m/z for PPM) ---
  tol_condition <- if (!is.null(tolerance_ppm)) {
    glue::glue("ABS(n.mz - t.mz) <= ({tolerance_ppm} / 1e6) * t.mz")
  } else {
    glue::glue("ABS(n.mz - t.mz) <= {tolerance_da}")
  }

  # VALUES list for DuckDB
  values_sql <- paste0("(", formatC(target_mz, digits = 12, format = "fg"), ")", collapse = ",")

  # --- query ---
  query <- glue::glue("
    WITH max_intensity AS (
      SELECT file, MAX(intensity) AS max_int
      FROM spectra
      {where_clause}
      GROUP BY file
    ),
    normalized AS (
      SELECT
        s.file,
        s.polarity,
        s.analyzer_type,
        s.ionization_type,
        s.collision_energy_level,
        s.collision_energy_voltage,
        s.mz,
        s.intensity,
        (s.intensity / m.max_int) * 100.0 AS rel_intensity
      FROM spectra s
      JOIN max_intensity m ON s.file = m.file
      {where_clause}
    ),
    matched AS (
      SELECT
        n.file,
        n.polarity,
        n.analyzer_type,
        n.ionization_type,
        n.collision_energy_level,
        n.collision_energy_voltage,
        COUNT(DISTINCT t.mz) AS match_count
      FROM normalized n
      JOIN (VALUES {values_sql}) AS t(mz)
        ON {tol_condition}
      WHERE n.rel_intensity >= {min_rel_intensity}
      GROUP BY n.file, n.polarity, n.analyzer_type, n.ionization_type,
               n.collision_energy_level, n.collision_energy_voltage
    )
    SELECT
      m.file,
      m.polarity,
      m.analyzer_type,
      m.ionization_type,
      m.collision_energy_level,
      m.collision_energy_voltage,
      m.match_count
    FROM matched m
    WHERE m.match_count = {length(target_mz)}
  ")

  res <- DBI::dbGetQuery(con, query)

  res$hmdb_id <- stringr::str_extract(res$file, "HMDB\\d+")

  res
}


find_peak_diff <- function(
    con, target_diff, tolerance = 0.01,
    polarity = NULL,
    collision_energy_level = NULL,
    voltage_range = NULL,
    min_rel_intensity = 5,
    mass_map = hmdb_mass_map
) {
  # 1. Create temporary table with hmdb_mass_map
  DBI::dbExecute(con, "DROP TABLE IF EXISTS hmdb_mass_map_tmp")
  DBI::dbWriteTable(con, "hmdb_mass_map_tmp", as.data.frame(mass_map), temporary = TRUE)

  # 2. Build optional filters
  meta_filters <- c("1=1")

  if (!is.null(polarity)) {
    meta_filters <- c(meta_filters, glue::glue("s1.polarity = '{polarity}'"))
  }
  if (!is.null(collision_energy_level)) {
    meta_filters <- c(meta_filters, glue::glue("s1.collision_energy_level = '{collision_energy_level}'"))
  }
  if (!is.null(voltage_range) && length(voltage_range) == 2) {
    meta_filters <- c(meta_filters,
                      glue::glue("s1.collision_energy_voltage BETWEEN {voltage_range[1]} AND {voltage_range[2]}"))
  }

  where_meta <- paste(meta_filters, collapse = " AND ")

  # 3. SQL query — now with more metadata
  query <- glue::glue("
    WITH max_intensity AS (
      SELECT file, MAX(intensity) AS max_int
      FROM spectra
      GROUP BY file
    ),
    normalized AS (
      SELECT
        s.file,
        s.polarity,
        s.analyzer_type,
        s.ionization_type,
        s.collision_energy_level,
        s.collision_energy_voltage,
        s.adduct_mass,
        s.mz,
        s.intensity,
        (s.intensity / m.max_int) * 100.0 AS rel_intensity
      FROM spectra s
      JOIN max_intensity m ON s.file = m.file
    ),
    precursors AS (
      SELECT
        f.file,
        f.polarity,
        f.analyzer_type,
        f.ionization_type,
        f.collision_energy_level,
        f.collision_energy_voltage,
        COALESCE(f.adduct_mass, m.adduct_mass) AS adduct_mass,
        CASE
          WHEN MIN(
            CASE
              WHEN p.rel_intensity >= {min_rel_intensity} THEN ABS(p.mz - COALESCE(f.adduct_mass, m.adduct_mass))
              ELSE NULL
            END
          ) <= {tolerance}
          THEN (
            SELECT p2.mz
            FROM normalized p2
            WHERE p2.file = f.file
              AND p2.rel_intensity >= {min_rel_intensity}
            ORDER BY ABS(p2.mz - COALESCE(f.adduct_mass, m.adduct_mass)) ASC
            LIMIT 1
          )
          ELSE COALESCE(f.adduct_mass, m.adduct_mass)
        END AS precursor_mz,
        CASE
          WHEN EXISTS (
            SELECT 1
            FROM normalized p
            WHERE p.file = f.file
              AND ABS(p.mz - COALESCE(f.adduct_mass, m.adduct_mass)) <= {tolerance}
              AND p.rel_intensity >= {min_rel_intensity}
          ) THEN 'real'
          ELSE 'virtual'
        END AS precursor_type
      FROM (
        SELECT DISTINCT
          s1.file,
          s1.polarity,
          s1.analyzer_type,
          s1.ionization_type,
          s1.collision_energy_level,
          s1.collision_energy_voltage,
          s1.adduct_mass,
          REGEXP_REPLACE(s1.file, '^.*(HMDB[0-9]+).*$', '\\1') AS hmdb_id
        FROM spectra s1
        WHERE {where_meta}
      ) f
      LEFT JOIN hmdb_mass_map_tmp m ON f.hmdb_id = m.database_id
      LEFT JOIN normalized p ON p.file = f.file
      GROUP BY f.file, f.polarity, f.analyzer_type, f.ionization_type,
               f.collision_energy_level, f.collision_energy_voltage,
               f.adduct_mass, m.adduct_mass
    )
    SELECT
      p.file,
      REGEXP_REPLACE(p.file, '^.*(HMDB[0-9]+).*$', '\\1') AS hmdb_id,
      p.precursor_mz,
      p.precursor_type,
      p.analyzer_type,
      p.ionization_type,
      p.collision_energy_level,
      p.collision_energy_voltage,
      f.mz AS frag_mz,
      f.intensity AS frag_intensity,
      (f.intensity / mi.max_int) * 100.0 AS frag_rel_intensity,
      ABS(p.precursor_mz - f.mz) AS diff_from_precursor
    FROM precursors p
    JOIN normalized f ON f.file = p.file
    JOIN max_intensity mi ON f.file = mi.file
    WHERE f.rel_intensity >= {min_rel_intensity}
      AND f.mz < p.precursor_mz
      AND ABS((p.precursor_mz - f.mz) - {target_diff}) <= {tolerance}
    ORDER BY p.file, diff_from_precursor
  ")

  df <- DBI::dbGetQuery(con, query)
  if ("precursor_type" %in% names(df)) {
    df$precursor_type <- factor(df$precursor_type, levels = c("real", "virtual"))
  }
  df
}


get_compound_name <- function(hmdb_id, hmdb_map) {
  hmdb_map$database_id <- trimws(hmdb_map$database_id)
  hmdb_id <- trimws(hmdb_id)
  idx <- match(hmdb_id, hmdb_map$database_id)
  name <- hmdb_map$compound_name[idx]
  ifelse(is.na(name) | name == "", hmdb_id, name)
}

get_compound_smiles <- function(hmdb_id, hmdb_map) {
  hmdb_map$database_id <- trimws(hmdb_map$database_id)
  hmdb_id <- trimws(hmdb_id)
  idx <- match(hmdb_id, hmdb_map$database_id)
  smiles <- hmdb_map$smiles[idx]
  ifelse(is.na(smiles) | smiles == "", hmdb_id, smiles)
}

smiles_to_svg <- function(smiles, width = 300, height = 300) {
  tf <- tempfile(fileext = ".svg")

  # Convert SMILES to SDFset
  molset <- ChemmineR::smiles2sdf(smiles)
  if (length(molset) == 0) return("")  # empty on failure

  # Extract the first molecule (SDF)
  mol <- molset[[1]]

  svglite::svglite(tf, width = width/64, height = height/64)
  on.exit(dev.off(), add = TRUE)
  ChemmineR::plot(mol)
  dev.off()

  paste(readLines(tf), collapse = "\n")
}

smiles_to_img_tag <- function(smiles, name, width = 150, height = 150) {
  # Directory inside www
  temp_dir <- file.path(app_sys("app/www/temp"))
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)

  # Safe filename
  file_name <- paste0(gsub("[^A-Za-z0-9]", "_", name), ".svg")
  file_path <- file.path(temp_dir, file_name)

  # Generate image if missing
  if (!file.exists(file_path)) {
    molset <- ChemmineR::smiles2sdf(smiles)
    if (length(molset) == 0) return("")
    mol <- molset[[1]]

    svglite::svglite(file_path, width = width / 40, height = height / 40)
    par(mar = c(3, 3, 3, 3))
    ChemmineR::plot(mol)
    dev.off()
  }

  # Return relative path so Shiny can serve it
  paste0("<img src='www/temp/", file_name, "' width='", width, "' height='", height, "'>")
}
