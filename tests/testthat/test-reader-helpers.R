test_that(".parse_logical accepts logical input", {
  # readxl and data.table::fread type a TRUE/FALSE column as logical, so the
  # most obvious spelling in a molecule definition file arrives here as a
  # logical and must not fall through to NULL.
  expect_true(pksimhelper:::.parse_logical(TRUE))
  expect_false(pksimhelper:::.parse_logical(FALSE))
  expect_null(pksimhelper:::.parse_logical(NA))

  expect_true(pksimhelper:::.parse_logical("yes"))
  expect_false(pksimhelper:::.parse_logical("No"))
  expect_true(pksimhelper:::.parse_logical(1))
  expect_null(pksimhelper:::.parse_logical("maybe"))
})

test_that(".id.molecule.cols ignores an empty additional matcher", {
  df <- data.frame(1, 2)
  colnames(df) <- c("Time [h]", "Sim A|Caffeine [umol/l]")

  # no additional matcher at all
  expect_equal(pksimhelper:::.id.molecule.cols(df, "caffeine")$cols, 2)

  # NA additional matcher, as produced by an empty "More" cell
  expect_equal(pksimhelper:::.id.molecule.cols(df, "caffeine", NA)$cols, 2)

  # a real additional matcher still narrows the search
  expect_equal(pksimhelper:::.id.molecule.cols(df, "caffeine", "sim a")$cols, 2)
  expect_null(pksimhelper:::.id.molecule.cols(df, "caffeine", "sim b", silent = TRUE))
})

test_that(".read.all.sheet returns one entry per CSV file", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  utils::write.csv(data.frame(`Time [h]` = c(0, 1), Value = c(1, 2), check.names = FALSE),
                   f, row.names = FALSE)

  res <- pksimhelper:::.read.all.sheet(basename(f), folder = dirname(f))
  expect_length(res, 1)
  expect_s3_class(res[[1]], "data.frame")
  expect_equal(colnames(res[[1]]), c("Time [h]", "Value"))
})
