################### is_single_string ###################

test_that("is_single_string works for single strings", {

  expect_true(is_single_string(""))
  expect_true(is_single_string(" "))
  expect_true(is_single_string("Hallo"))
  expect_true(is_single_string("1234"))
})

test_that("is_single_string works for vectors", {

  expect_true(is_single_string(c("")))
  expect_false(is_single_string(c("A", "B")))
})

test_that("is_single_string works for lists", {

  expect_false(is_single_string(list("")))
  expect_false(is_single_string(list(" ", " ")))
})

test_that("is_single_string works for NA and NULL", {

  expect_false(is_single_string(NULL))
  expect_false(is_single_string(NA))
  expect_false(is_single_string(NA_character_))
  expect_false(is_single_string(NA_complex_))
  expect_false(is_single_string(NA_real_))
})

test_that("is_single_string works for not-empty-allowed", {

  expect_false(is_single_string("", can.be.empty = F))

  expect_true(is_single_string(" ", can.be.empty = F, trim.check = F))
  expect_false(is_single_string(" ", can.be.empty = F, trim.check = T))
})

################### is_single_logical ###################

test_that("is_single_logical works for one logical input", {

  expect_true(is_single_logical(T))
  expect_true(is_single_logical(F))
})

test_that("is_single_logical works for one non-logical input", {

  expect_false(is_single_logical(""))
  expect_false(is_single_logical(1))
  expect_false(is_single_logical(0))
  expect_false(is_single_logical(1.2))
  expect_false(is_single_logical(c()))
  expect_false(is_single_logical(list()))
  expect_false(is_single_logical(data.frame()))
})

test_that("is_single_logical works for NULL and NA", {

  expect_false(is_single_logical(NA))
  expect_false(is_single_logical(NULL))
  expect_false(is_single_logical(NA_character_))
  expect_false(is_single_logical(NA_complex_))
  expect_false(is_single_logical(NA_integer_))
})

test_that("is_single_logical works for vectors and lists", {

  expect_false(is_single_logical(c(T, F)))
  expect_false(is_single_logical(list(T, F)))
  expect_false(is_single_logical(list(T)))
  expect_true(is_single_logical(c(T)))
})

################### is_single_numeric ###################

test_that("is_single_numeric works for one numeric input", {

  expect_true(is_single_numeric(1))
  expect_true(is_single_numeric(1.2))
})

test_that("is_single_logical works for one non-numeric input", {

  expect_false(is_single_numeric(""))
  expect_false(is_single_numeric(F))
  expect_false(is_single_numeric(T))
  expect_false(is_single_numeric(c()))
  expect_false(is_single_numeric(list()))
  expect_false(is_single_numeric(data.frame()))
})

test_that("is_single_logical works for NULL and NA", {

  expect_false(is_single_numeric(NA))
  expect_false(is_single_numeric(NULL))
  expect_false(is_single_numeric(NA_character_))
  expect_false(is_single_numeric(NA_complex_))
  expect_false(is_single_numeric(NA_integer_))
  expect_false(is_single_numeric(NA_real_))
})

test_that("is_single_numeric works for vectors and lists", {

  expect_false(is_single_numeric(c(1, 2)))
  expect_false(is_single_numeric(list(1, 2)))
  expect_false(is_single_numeric(list(1)))
  expect_true(is_single_numeric(c(1)))
})

################### has_units ###################

test_that("has_units works for units", {

  a <- units::as_units(1, "m")
  c <- units::as_units(c(1, 5), "m")
  b <- list(a = c, b = c)
  d <- data.frame(a = c, b = c)

  expect_true(has_units(a))
  expect_true(has_units(b))
  expect_true(has_units(c))
  expect_true(has_units(d))
})

test_that("has_units works for non-units", {

  a <- 1
  c <- c(1, 5)
  b <- list(a = c, b = c)
  d <- data.frame(a = c, b = c)

  expect_false(has_units(a))
  expect_false(has_units(b))
  expect_false(has_units(c))
  expect_false(has_units(d))
})


test_that("has_units works for NA and NULL", {

  expect_false(has_units(NA))
  expect_false(has_units(NULL))
})

################### is_string_list ###################

test_that("is_string_list works for valid input", {

  expect_true(is_string_list("A"))
  expect_true(is_string_list(c("A", "B")))
  expect_true(is_string_list(list("A", "B")))
  expect_true(is_string_list(list("A")))
})

test_that("is_string_list works for string like", {

  expect_true(is_string_list(1, string.like = T))
  expect_true(is_string_list(c(12.4, 43), string.like = T))
  expect_true(is_string_list(list(12.4, T), string.like = T))
})

test_that("is_string_list works for empty strings", {

  expect_true(is_string_list("", can.be.empty = T))
  expect_true(is_string_list(c(" ", ""), can.be.empty = T))
})

test_that("is_string_list works for invalid input", {

  expect_false(is_string_list(""))
  expect_false(is_string_list(c("A", " ")))
  expect_false(is_string_list(list("A", NA)))
  expect_false(is_string_list(list("A", NA)))
  expect_false(is_string_list(list(1.2, 5)))
})

################### is_unit ###################

test_that("is_unit works", {

  expect_true(is_unit("m"))
  expect_true(is_unit("m/s"))
  expect_false(is_unit("Foo"))
  expect_false(is_unit(""))
  expect_false(is_unit(" "))
  expect_false(is_unit(list()))
  expect_false(is_unit(c("m", "sec")))
})

################### is_dose_unit ###################

test_that("is_dose_unit works", {

  expect_true(is_dose_unit("kg"))
  expect_true(is_dose_unit("µg"))
  expect_true(is_dose_unit("mg"))
  expect_true(is_dose_unit("ng"))
  expect_true(is_dose_unit("g"))
  expect_true(is_dose_unit("pg"))
  expect_true(is_dose_unit("kg"))
  expect_false(is_dose_unit("m/s"))
  expect_false(is_dose_unit("Foo"))
  expect_false(is_dose_unit(""))
  expect_false(is_dose_unit(" "))
  expect_false(is_dose_unit(list()))
  expect_false(is_dose_unit(c("m", "sec")))
})


################### has_dose ###################

test_that("has_dose works", {

  expect_true(has_dose(units::as_units(12, "mg"))) # true
  expect_true(has_dose(units::as_units(12.4, "µg"))) # true
  expect_false(has_dose(units::as_units(12.4, "m"))) # false
  expect_false(has_dose(12.4)) # false
})


################### to_range ###################

test_that("to_range works", {

  expect_equal(to_range(NA_real_, NA_real_, NA), c(NA_real_, NA_real_)) # will return c(NA, NA)
  expect_equal(to_range(NA, 1, "mg"), units::as_units(c(NA, 1), "mg"))
  expect_equal(to_range(NA, "1", "mg"), units::as_units(c(NA, 1), "mg"))
  expect_error(to_range(1, NA, NA)) # error
})
