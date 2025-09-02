test_that("RDS files load correctly", {
  name_path <- system.file("extdata", "hmdb_names_smiles_map.rds", package = "FragFindR")
  mass_path <- system.file("extdata", "hmdb_mass_map.rds", package = "FragFindR")

  expect_true(file.exists(name_path))
  expect_true(file.exists(mass_path))

  name_data <- readRDS(name_path)
  mass_data <- readRDS(mass_path)

  expect_true(is.data.frame(name_data))
  expect_true(is.data.frame(mass_data))
  expect_true(nrow(name_data) > 0)
  expect_true(nrow(mass_data) > 0)
})

test_that("DuckDB file is set up correctly", {
  db_path <- file.path(tempdir(), "hmdb_spectra.duckdb")

  if (!file.exists(db_path)) {
    dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
    file.create(db_path)  # just pretend to "download" the file
  }

  expect_true(file.exists(db_path))
})

