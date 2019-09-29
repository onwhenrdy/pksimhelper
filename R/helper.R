
####################################################
# Input Helper
####################################################

#' Convinience Helper to Test for any `NA` value
#'
#' @param x The test value.
#'
#' @return `TRUE` if any value in `x` is `NA`.
#' @export
#'
#' @family helper functions
#'
any_na <- function(x) {
  any(is.na(x))
}

#' Convinience Helper to Test for single `NA` value
#'
#' @param x The test value.
#'
#' @return `TRUE` if first value in `x` is `NA` and `length(x) == 1`.
#' @export
#'
#' @family helper functions
#'
single_na <- function(x) {
  any(is.na(x)) && length(x) == 1
}

#' Convinience Helper to Test for any `NULL` value
#'
#' @param x The test value.
#'
#' @return `TRUE` if any value in `x` is `NULL`.
#' @export
#'
#' @family helper functions
#'
any_null <- function(x) {
  any(is.null(x))
}


#' Checks if Input is a Single String
#'
#' @param input The input that should be tested
#' @param can.be.empty If the input (if string) is allowed to be empty.
#' @param trim.check If the input (if string) should be trimmed before the check.
#'
#' @return `TRUE` if the input is a single string of vector with a single string.
#' Lists with one string element will return `FALSE`.
#'
#' @family helper functions
#'
#' @export
#'
is_single_string <- function(input, can.be.empty = T, trim.check = T) {

  is.single.str <- (is.character(input) && length(input) == 1)

  if (trim.check)
    input <- base::trimws(input)

  return(is.single.str && (can.be.empty || nchar(input) > 0 ) && !is.na(input))
}


#' Tests if Input is a String or List-like Object of Strings
#'
#' @param vec The input object.
#' @param string.like If input is allowed to be convertible to string for the test.
#' @param can.be.empty If empty strings are allowed.
#'
#' Returns `TRUE` for empty lists. For lists with NA-values this function will return `FALSE`.
#'
#' @return `TRUE`, if the input is a string or list-like object of strings, else `FALSE`.
#' @export
#'
#' @family helper functions
#'
is_string_list <- function(vec, string.like = F, can.be.empty = F) {

  if (any(is.na(vec)))
    return(F)

  if (string.like)
    vec <- paste(vec)

  all(sapply(vec, is_single_string, can.be.empty = can.be.empty))
}


#' Checks if Input is a Single Logical Value
#'
#' @param input The input that should be tested
#'
#' @return `TRUE` if the input is a single logical value. Lists with one logical value will return `FALSE`.
#' @export
#'
#' @family helper functions
#'
is_single_logical <- function(input) {

  is.logical(input) && length(input) == 1 && !is.na(input)
}


#' Checks if Input is a Single Bumeric Value
#'
#' @param input The input that should be tested
#'
#' @return `TRUE` if the input is a single numeric value. Lists with one numeric value will return `FALSE`.
#' @export
#'
#' @family helper functions
#'
is_single_numeric <- function(input) {

  is.numeric(input) && length(input) == 1 && !is.na(input)
}

# TODO: Tests
#' Tests if Input is a List-like Object of Numerics
#'
#' @param vec The input object.
#' @param can.be.empty If NA-values are allowed.
#'
#' Returns `TRUE` for empty lists.
#'
#' @return `TRUE`, if the input is a list-like object of numerics, else `FALSE`.
#' @export
#'
#' @family helper functions
#'
is_numeric_list <- function(vec, can.be.empty = F) {

  if (length(vec) == 0)
    return(F)

  if (!can.be.empty && any(is.na(vec)))
    return(F)

  all(sapply(vec, is_single_numeric))
}



#' Tests if a Data Structure has Units Attached.
#'
#' @param input The data structure that should be tested.
#'
#' @return `TRUE` if the data structure inherits from `units` (must use the `units` package).
#' @export
#'
#' @family helper functions
#'
has_units <- function(input) {

  if (length(input) == 0)
    return(FALSE)

  test_fn <- function(x) inherits(x, "units")
  if (test_fn(input))
    return(TRUE)

  all(sapply(input, test_fn))
}

#' Tests for a Valid Dose Unit String
#'
#' This function will only work correctly for lower case strings (e.g. kg and not KG).
#' For non-strings, NA or NULL this function will return False.
#'
#' @param str The string that should be tested.
#'
#' @return `TRUE`, if a valid dose unit is provided, else `FALSE`.
#' @export
#'
#' @family helper functions
#'
#' @examples
#'
#' is_dose_unit("") # F
#' is_dose_unit(NA) # F
#' is_dose_unit("mg") # T
#' is_dose_unit("pg") # T
#' is_dose_unit("µg") # T
#' is_dose_unit("µg/kg") # F
#'
is_dose_unit <- function(str) {

  if (!is_single_string(str, can.be.empty = F)) {
    return(F)
  }

  test <- units::as_units(1, "mg")
  tryCatch({units(test) <- str; T},
           error = function(e) F,
           warning = function(w) F)
}

#' Tests for a Valid Time Unit String
#'
#' This function will only work correctly for lower case strings (e.g. h and not H).
#' For non-strings, NA or NULL this function will return False.
#'
#' @param str The string that should be tested.
#'
#' @return `TRUE`, if a valid time unit is provided, else `FALSE`.
#' @export
#'
#' @family helper functions
#'
#' @examples
#'
#' is_time_unit("") # F
#' is_time_unit(NA) # F
#' is_time_unit("h") # T
#' is_time_unit("sec") # T
#' is_time_unit("weeks") # T
#'
is_time_unit <- function(str) {

  if (!is_single_string(str, can.be.empty = F)) {
    return(F)
  }

  test <- units::as_units(1, "h")
  tryCatch({units(test) <- str; T},
           error = function(e) F,
           warning = function(w) F)
}

#' Tests if a Single Numeric Value has a valid Dosing Unit attached
#'
#' @param value A single numeric value.
#'
#' @return `TRUE` if it is a single numeric value with a valid dosing unit (e.g. mg, µg, ...)
#'   attached, else `FALSE`.
#' @export
#'
#' @family helper functions
#'
#' @examples
#'
#' has_dose(units::as_units(12, "mg")) # true
#' has_dose(units::as_units(12.4, "\U00B5g")) # true
#' has_dose(units::as_units(12.4, "m")) # false
#' has_dose(12.4) # false
#'
has_dose <- function(value) {

  if (!is_single_numeric(value) || !has_units(value)) {
    return(F)
  }


  tryCatch({units(value) <- "mg"; T},
           error = function(e) F,
           warning = function(w) F)
}


#' Checks if a String is a Valid Unit
#'
#' @param str The input string.
#'
#' @return `TRUE` if a valid unit, else `FALSE`.
#' @export
#'
#' @family helper functions
#'
is_unit <- function(str) {

  if (!is_single_string(str, can.be.empty = F))
    return(F)

  ok <- tryCatch({
    units::as_units(str)
    T
  },
  error = function(e) {
    F
  })
  return(ok)
}


#' Converts two Values and a Unit String to a Valid Range
#'
#' @param min The minumum value (will be converted by `as.numeric`).
#' @param max The maximum value (will be converted by `as.numeric`)
#' @param unit.str A unit string (e.g. mg).
#'
#' @return A vector with `length` 2 that has a unit attached. If `min` and `max` are both `NA` no unit
#' will be attached.
#'
#' If `min` and `max` are both `NA` `unit.str` is ignored and will not be evaluated.
#'
#' @family helper functions
#'
#' @export
#'
#' @examples
#'
#' to_range(NA, NA, NA) # will return c(NA, NA)
#' to_range(NA, 1, "mg") # will return c(NA, 1) with attached unit "mg"
#' to_range(NA, "1", "mg") # will return c(NA, 1) with attached unit "mg"
#' # to_range(1, NA, NA) # error
#'
to_range <- function(min, max, unit.str) {

  result <- as.numeric(c(min, max))
  if (purrr::every(result, is.na))
    return(result)

  if (!is_unit(unit.str))
    stop(paste("unit string <", unit.str, "> is not a valid unit string"), call. = F)

  units::as_units(result, unit.str)
}


# TODO: Tests and docs
#' Converts Values from one unit into another.
#'
#' @param values Input values that should be converted from `from` to `to`. Input can be a vector/`tibble`
#'   or `data.frame`.
#' @param from  Unit string that is compatible with `to`.
#' @param to f Unit string that is compatible with `from`.
#' @param attach.unit If `TRUE` the unit `from` is attached to the returned values.

#' If values already has a unit attached `from` is ignored.
#'
#' @family helper functions
#'
#' @return A
#' @export
#'
#' @examples
#'
convert_values <- function(values, from = NA, to, attach.unit = F) {

  if (!is_unit(to))
    stop("`to` is not a valid unit.", call. = F)

  if (!has_units(values)) {
    if (!is_unit(from))
      stop("`from` is not a valid unit.", call. = F)

    units(values) <- from
  }

  units(values) <- to

  if (!attach.unit)
    values <- rm_units(values)

  return(values)
}

#' Removes Units from an Object
#'
#' If the input does not have units attached it is returned.
#'
#' @param x The input object.
#'
#'
#' @return `x` without units attached.
#'
#' @family helper functions
#'
#' @export
#'
rm_units <- function(x) {
  if (has_units(x))
    x <- units::drop_units(x)

  return(x)
}










# checks if input is a single character vector
.is.single.string <- function(input) {

  is.character(input) && length(input) == 1
}

####################################################
# FILE HELPER
####################################################

.read.csv <- function(file, sep, dec, header,
                      encoding = "UTF-8",
                      strip.white = TRUE) {

  as.data.frame(data.table::fread(file = file,
                    header = header,
                    sep = sep,
                    dec = dec,
                    encoding = encoding,
                    strip.white = strip.white))
}


# gather ids from a profile list (e.g. observed data list)
.gather.ids <- function(profile.list, ids) {
  results = list()

  for (profile in profile.list) {
    if (!is.profile(profile))
      stop("Entry in profile.list is not a profile")

    if (profile$id %in% ids)
      results <- append(results, list(profile))
  }

  return(results)
}

# checks if file exists
.does.file.exist <- function(file) {

  if (!.is.single.string(file))
    stop(paste("File <", file, "> must be a character vector"))

  if (file.access(file, mode = 0) == 0) TRUE else FALSE
}

# checks if a file is readble
.is.file.readable <- function(file) {

  if (!.is.single.string(file))
    stop(paste("File <", file, "> must be a character vector"))
  #TODO
  return(TRUE)
  if (file.access(file, mode = 4) == 0) TRUE else FALSE
}

# identifies a file extension
# returns tsv, csv, xls ord xlsx
.identify.file.ext <- function(file) {
  if (!.does.file.exist(file))
    stop(paste("file <", file, "> does not exist"))

  ext <- tolower(tools::file_ext(file))
  if (ext == "csv" || ext == "tsv")
    return(ext)

  if (ext == "xls" || ext == "xlsx")
    return("xls")

  return(NA)
}

# returns a list of sheets from an excel file
.sheets <- function(file) {
  return(readxl::excel_sheets(file))
}

# checks if an excel sheet exists in file
.does.sheet.exist <- function(file, sheet) {

  sheets <- .sheets(file)
  if (is.numeric(sheet) && length(sheets) >= sheet)
    return(TRUE)
  else
    return(is.element(sheet, sheets))

  return(FALSE)
}

####################################################
# UNIT HELPER
####################################################

# returns a unit from a string (e.g. from "Time [h]")
.extract.unit <- function(str) {

  pattern <- "(\\[.*?\\])"
  matches <- gregexpr(pattern, str)
  overlap <- regmatches(str, matches)
  overlap_clean <- unlist(overlap)
  overlap_clean <- sub(".*\\[(.*)\\].*", "\\1", overlap_clean, perl = TRUE)


  if (identical(overlap_clean, character(0)))
    stop(paste("Could not find unit in string <", str , ">"))

  if (tolower(overlap_clean) == "bpm")
    return(units::as_units("beats/min"))

  unit <- tryCatch({
    units::as_units(overlap_clean)
  },
  error = function(e) {
    stop(paste("String <", str , "> is not a valid unit"))
  })

  return(unit)
}

.has.mol.unit <- function(unit) {
  grepl("mol", tolower(units::deparse_unit(units::as_units(units(unit)$numerator))), fixed = TRUE)
}

.has.mass.unit <- function(unit) {
  grepl("g", tolower(units::deparse_unit(units::as_units(units(unit)$numerator))), fixed = TRUE)
}

# converts data vector or single number with unit from to unit to
.convert.units <- function(data, from.unit, to.unit, MW = NA) {

  if (is.data.frame(data))
    data <- as.matrix(data)

  units(data) <- from.unit

  any.mass <- .has.mass.unit(from.unit) || .has.mass.unit(to.unit)
  if (any.mass) {
    from.mol <- .has.mol.unit(from.unit)
    to.mol <- .has.mol.unit(to.unit)

    if ((xor(from.mol, to.mol))) {
      # need special conversion
      if (is.na(MW))
        stop("Conversion needs molecular weight")

      if (from.mol) {
        # -> to ?g/?
        units(data) <- from.unit # ?mol/?
        units(data)$numerator <- "mol" # mol/?
        data <- data * MW

        target <- units::as_units("g")
        denom <- units(data)$denominator
        if (length(denom) > 0)
          target <- target / units::as_units(denom)
        data <- units::drop_units(data)
        units(data) <- target
      } else {
        # -> to ?mol/?
        units(data) <- from.unit # ?g/?
        units(data)$numerator <- "g" # g/?
        data <- data / MW

        target <- units::as_units("mol")
        denom <- units(data)$denominator
        if (length(denom) > 0)
          target <- target / units::as_units(denom)
        data <- units::drop_units(data)
        units(data) <- target
      }
    }
  }

  units(data) <- to.unit

  if (is.matrix(data))
    data <- as.data.frame(data)

  return(units::drop_units(data))
}

####################################################
# MISC HELPER
####################################################

# adds alpha to a colour
add.alpha <- function(cols, alpha) grDevices::rgb(t(grDevices::col2rgb(cols)/255), alpha = alpha)

# finds a molecule with id in list of molecules
.find.molecule.from.id <- function(molecules, id) {

  match.id <- tolower(trimws(id))
  for (mol in molecules) {
    mol.id <- tolower(trimws(mol$id))
    if (mol.id == match.id)
      return(mol)
  }

  return(NA)
}

is.valid <- function(x,...) {
  if (is.profile(x))
    return(is.valid.profile(x, ...))
  else if (is.matched.profiles(x)) {
    result <- unlist(lapply(x$profiles, is.valid.profile, ...))
    return(all(result))
  }else if (is.list(x)) {
    result <- unlist(lapply(x, is.valid, ...))
    return(all(result))
  }
}

`%notin%` <- Negate(`%in%`)

is.valid.profile <- function(profile, msg = c("warning", "message", "error", "noop")) {

  msg.type <- match.arg(msg)
  report <- function(x) {}
  if (msg.type == "warning")
    report = warning
  else if (msg.type == "message")
    report = message
  else if (msg.type == "error")
    report = stop

  error = FALSE
  id <- profile$id

  if (!is.profile(profile))
    stop("Input must be of class profile")

  if (!is.molecule(profile$molecule)) {
    error = TRUE
    report(paste("Found profile with no molecule (id: ", id, ")"))
  }

  if (!is.character(profile$reference)) {
    error = TRUE
    report(paste("Found profile with no character reference (id: ", id, ")"))
  }

  types <- c("individual", "population")
  d.types <- c("mean", "individual")
  origins <- c("sim", "obs")
  if (profile$origin %notin% origins) {
    error = TRUE
    report(paste("Found profile with unknown origin (id: ", id, ")"))
  }

  if (!(profile$type %in% types)) {
    error = TRUE
    report(paste("Found profile with unknown type (id: ", id, ")"))
  }

  if (!(profile$data.type %in% d.types)) {
    error = TRUE
    report(paste("Found profile with unknown data.type (id: ", id, ")"))
  }

  tryCatch({test <- units(profile$time.unit)}, error = function(e) {
    error = TRUE
    report(paste("Found profile with no time unit (id: ", id, ")"))
  })

  tryCatch({test <- units(profile$value.unit)}, error = function(e) {
    error = TRUE
    report(paste("Found profile with no value unit (id: ", id, ")"))
  })

  data <- profile$data
  names <- colnames(data)
  if ("Time" %notin% names) {
    error = TRUE
    report(paste("Found profile data with no time column (id: ", id, ")"))
  }

  if (is.unsorted(data$Time, strictly = TRUE)) {
    error = TRUE
    report(paste("Found profile data with unsorted time column (id: ", id, ")"))
  }

  if (profile$data.type == "mean") {
    if ("Avg" %notin% names) {
      error = TRUE
      report(paste("Found mean profile data with no Avg column (id: ", id, ")"))
    }

    if ("Min" %notin% names) {
      error = TRUE
      report(paste("Found mean profile data with no Min column (id: ", id, ")"))
    }

    if ("Max" %notin% names) {
      error = TRUE
      report(paste("Found mean profile data with no Nax column (id: ", id, ")"))
    }

    if (any(is.nan(data$Avg)) || any(is.na(data$Avg))) {
      error = TRUE
      report(paste("Found mean profile data NA/NaN Avg data (id: ", id, ")"))
    }

    min.na <- sum(is.na(data$Min)) + sum(is.nan(data$Min))
    max.na <- sum(is.na(data$Max)) + sum(is.nan(data$Max))

    if (min.na != max.na) {
      error = TRUE
      report(paste("Found unbalanced NA/NaN data in Min/Max columns (id: ", id, ")"))
    }

    if (min.na != 0 && min.na != length(data$Min)) {
      error = TRUE
      report(paste("Found mixed NA/NaN data in Avg/Min/Max columns (id: ", id, ")"))
    }
  } else  {
    # individual -> no NA/NaNs allowed
    if (sum(apply(data, 2, is.nan)) >=1 || sum(apply(data, 2, is.na)) >=1) {
      error = TRUE
      report(paste("Found NA/NaN data for indivdual profile (id: ", id, ")"))
    }
  }
  return(!error)
}

.has.units <- function(obj) {
  inherits(obj, "units")
}

.rm.units <- function(obs) {
  if (.has.units(obs))
    obs <- units::drop_units(obs)

  return(obs)
}
