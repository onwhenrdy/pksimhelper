#' Pipe
#'
#' Put description here
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs specify what lhs and rhs are
#' @examples
#' # some examples if you want to highlight the usage in the package
NULL


.ls_override <- function(default, new) {

  res <- c(new, default)
  res[!base::duplicated(names(res))]
}


# checks if input is a single character vector
.is.single.string <- function(input) {

  is.character(input) && length(input) == 1
}

is.blank <- function(x){
  if (is.function(x)) return(FALSE)

  return(
    is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)))
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
  if (length(ids) == 1 && is.na(ids))
    return(list())

  results = sapply(ids,function(x) NULL)

  for (profile in profile.list) {
    if (!is.profile(profile))
      stop("Entry in profile.list is not a profile")

    if (profile$id %in% ids)
      results[[profile$id]] <- profile
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

  if (tolower(overlap_clean) == "day(s)")
    return(units::as_units("days"))

  unit <- tryCatch({
    units::as_units(overlap_clean)
  },
  error = function(e) {
    stop(paste("String <", str , "> is not a valid unit"))
  })

  return(unit)
}

.has.units <- function(obj) {
  inherits(obj, "units")
}

.rm.units <- function(obs) {
  if (.has.units(obs))
    obs <- units::drop_units(obs)

  return(obs)
}

.has.mol.unit <- function(unit) {
  if (units::deparse_unit(unit) == "")
    return(F)

  grepl("mol", tolower(units::deparse_unit(units::as_units(units(unit)$numerator))), fixed = TRUE)
}

.has.mass.unit <- function(unit) {
  if (units::deparse_unit(unit) == "")
    return(F)

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


# sanatizes files names
.sanatize_filename <- function(file) {

  file <- gsub(" +", "", file)
  file <- gsub("/", "_", file)
  file <- gsub("&", "+", file)
  iconv(file, "latin1", "ASCII", sub = "x")
}


# parses a string or number and tries to convert to logical (TRUE/FALASE)
# Returns NULL if not possible
.parse_logical <- function(input) {

  if (is.character(input)) {

    input <- trimws(tolower(input))

    if (input == "t" || input == "true" ||
        input == "yes" || input == "1")
      return(TRUE)

    if (input == "f" || input == "false" ||
        input == "no" || input == "0")
      return(FALSE)
  }
  else if (is.numeric(input)) {

    if (input == 1)
      return(TRUE)

    if (input == 0)
      return(FALSE)
  }

  return(NULL)
}



#' Tests if a matched profiles has molecules with fraction unit
#'
#' @param profile A matched profile
#'
#' @return TRUE if the molecule of any profile has fraction unit else FALSE
#'
has_fraction <- function(profile) {

  for (pro in profile$profiles) {

    if (pro$mol$is.fraction)
      return(TRUE)

  }

  return(FALSE)
}




