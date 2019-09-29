################### Parpser Helper ###################

#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @family parser functions
#'
#' @examples
file_exists <- function(file) {

  if (!is_single_string(file))
    stop(paste("File <", file, "> must be a character vector"), call. = FALSE)

  if (file.access(file, mode = 0) == 0) TRUE else FALSE
}


#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @family parser functions
#'
#' @examples
file_type <- function(file) {
  if (!file_exists(file))
    stop(paste("file <", file, "> does not exist"), call. = FALSE)

  ext <- tolower(tools::file_ext(file))
  if (ext == "csv" || ext == "tsv")
    return(ext)

  if (ext == "xls" || ext == "xlsx")
    return("excel")

  return(NA)
}


#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @family parser functions
#'
#' @examples
excel_sheets <- function(file) {
  if (!file_exists(file))
    stop(paste("file <", file, "> does not exist"), call. = FALSE)

  return(readxl::excel_sheets(file))
}

#' Title
#'
#' @param file
#' @param file.type
#' @param excel.sheet
#' @param csv.sep
#' @param csv.encoding
#' @param csv.dec
#' @param csv.group
#' @param csv.comment
#' @param skip.lines
#'
#' @return
#' @export
#'
#' @family parser functions
#'
#' @examples
read_file_table <- function(file,
                            file.type = c("auto", "excel", "csv", "tsv"),
                            excel.sheet = 1,
                            csv.sep = if (file.type == "tsv") "\t" else ",",
                            csv.encoding = "guess",
                            csv.dec = ".",
                            csv.group = ",",
                            csv.comment = "",
                            skip.lines = 0) {

  file.type = match.arg(file.type)
  if (!file_exists(file))
    stop(paste("File <", file, "> does not exist"), call. = F)

  if (file.type == "auto")
    file.type <- file_type(file)

  if (is.na(file.type))
    stop(paste("Could not determine file type for file <", file, ">"), call. = F)

  # csv
  if (file.type == "csv" || file.type == "tsv") {

    locale = readr::default_locale()
    if (csv.encoding != "guess") {
      locale$encoding = csv.encoding
    }
    locale$decimal_mark = csv.dec
    locale$grouping_mark = csv.group

    return(suppressMessages(readr::read_delim(file,
                             comment = csv.comment,
                             delim = csv.sep,
                             locale = locale,
                             trim_ws = TRUE,
                             progress = FALSE,
                             skip_empty_rows = TRUE,
                             skip = skip.lines)))
  }

  suppressMessages(readxl::read_excel(file,
                                      skip = skip.lines,
                                      sheet = excel.sheet,
                                      col_names = TRUE))
}


#' Extract a Unit from an Input String
#'
#' @param str That string that should contain the unit.
#' @param extra.ident.fn Extra identification function.
#'
#' The input string `str` is searched for the pattern "[X]" with "X" being the unit.
#' If no pattern is found `NA` is returned. If the unit is unknown to `units` an error is raised.
#'
#' The user can provide a function `extra.ident.fn` with the signature `function(x) = ...`.
#' If the unit cannot be identified `extra.ident.fn` is invoked (if not `NULL`).
#' The function should return a unit as `units::as_units(...)` or `stop` if it cannot be identified.
#' The pattern inside `[]` will be trimmed and converted to lower case.
#'
#' @family parser functions
#'
#' @return `NA` or the unit as `units::as_units`.
#' @export
#'
#' @examples
#'
#' extract_unit("[DomInik]",
#'              function(x) if(x == "dominik") units::as_units("h") else stop("Oh noo"))
#'
#' extract_unit("Bla bla [h]") # returns as_units("h")
#'
extract_unit <- function(str, extra.ident.fn = NULL) {

  pattern <- "(\\[.*?\\])"
  matches <- gregexpr(pattern, str)
  overlap <- regmatches(str, matches)
  overlap_clean <- unlist(overlap)
  overlap_clean <- sub(".*\\[(.*)\\].*", "\\1", overlap_clean, perl = TRUE)


  if (identical(overlap_clean, character(0)))
    return(NA)

  overlap_clean <- tolower(trimws(overlap_clean))
  if (overlap_clean == "bpm")
    return(units::as_units("beats/min"))

  unit <- tryCatch({
    units::as_units(overlap_clean)
  },
  error = function(e) {

    if (!is.null(extra.ident.fn))
    {
      return(extra.ident.fn(overlap_clean))
    }
    stop(paste("String <", str , "> has not a valid unit"))
  },
  warning = function(w) {

    if (!is.null(extra.ident.fn))
    {
      return(extra.ident.fn(overlap_clean))
    }
    stop(paste("String <", str , "> has not a valid unit"))
  }
  )

  return(unit)
}


#' Identifies and Returns a Match in a `tibble`/`data.frame` Header or String Vector
#'
#' @param x `tibble`/`data.frame` or String Vector
#' @param match A single string that should be identified
#' @param multi.matches If `TRUE` multiple matches are returned and not only the best match.
#'
#' @family parser functions
#'
#' The function tries to match as follows:
#'   1. global match (has the match pattern)
#'   2. If ambiguous: Match strings that start with `match`
#'   3. If ambiguous: Match strings that start with `match` and have whitespace the match.
#'
#' If `multi.matches` is set to `TRUE` only option (1) will be applied.
#'
#' @return The position in the string vector or header or `NA` is not found or ambiguous.
#' @export
#'
ident_column <- function(x, match, multi.matches = F) {

  if (!is_single_string(match))
    stop("match must be a single non-empty string", call. = FALSE)

  cols <- x
  if (tibble::is_tibble(x) || is.data.frame(x))
    cols <- colnames(x)

  cols <- tolower(trimws(cols))

  # global match
  all <- which(stringr::str_detect(cols, match))
  if (length(all) == 1) {
    return(all[1])
  }

  # only begin with match
  begin <- which(stringr::str_detect(cols, paste0("^", match)))
  if (length(begin) == 1) {
    return(begin[1])
  }

  if (multi.matches && length(begin) >= 1)
    return(begin)

  # only begin with whitespace afterwards
  begin_ws <- which(stringr::str_detect(cols, paste0("^", match, "\\s")))
  # perfect match
  if (length(begin_ws) == 1) {
    return(begin_ws[1])
  }
  return(NA)
}


#' Tries to Identify a Time Entry in a `tibble`/`data.frame` Header or String Vector
#'
#' @param df A `tibble`/`data.frame` or string vector.
#' @param matches String or vector of strings that are matched one after another.
#'
#' The function tries to match all patterns defined in `matches` first.
#'   If no entry is found or ambiguous matches it searches for units ("[X]") that are convertible
#'   to a valid time unit.
#'
#' @return NA if no non-ambiguous entry could be identified, else the entry index.
#' @export
#'
#' @family parser functions
#'
#' @importFrom magrittr %>%
#'
ident_time_column <- function(df, matches = c("time", "tad", "t")) {

  # 1. search for matches
  # 2. search for [x] that is convertible to time unit

  # 1
  for (match in matches) {
    col <- ident_column(df, match)
    if (!is.na(col))
      return(col)
  }

  # 2
  cols <- df
  if (tibble::is_tibble(df) || is.data.frame(df))
    cols <- colnames(df)

  cols <- tolower(trimws(df))

  cols <- cols %>%
              purrr::map(extract_unit) %>%
              function(x) if (is.na(x)) FALSE else is_time_unit(x) %>%
              which
  if (length(cols) == 1)
    return(cols)

  return(NA)
}



################### Master Parser ###################

#' Title
#'
#' @param file_or_data
#' @param id
#' @param molecules
#' @param entry.sep
#' @param ...
#'
#' @return
#' @export
#'
#' @family parser functions
#'
#' @importFrom magrittr %>%
#'
#' @examples
#'
parse_master <- function(file_or_data,
                         id,
                         molecules,
                         entry.sep = ",",
                         ...) {

  if (tibble::is_tibble(file_or_data) || is.data.frame(file_or_data)) {
    data <- file_or_data
    ref <- "data.frame"
  } else {
    data <- read_file_table(file_or_data, ...)
    ref <- basename(file_or_data)
  }

  rep.fn <- function(entry, name) {
    if (!any_na(entry))
      message(paste("Identified", name, "in column", entry))
    else
      message(paste("WARNING: Could not identify", name, "column"))
  }

  # pop and sims
  ############################################################################
  # identify columns (pop id, pop name, population id, population name)
  pop_id <- ident_column(data, "pop.*(id|identifier|name)$")

  # identify columns (pop file, population file)
  pop_file <- ident_column(data, "pop.*(file)$")

  # identify columns (sim id, sim name, simulation id, simulation name)
  sim_id <- ident_column(data, "sim.*(id|identifier|name)$")

  # identify columns (sim file, simulation file)
  sim_file <- ident_column(data, "sim.*(file)$")

  pop_file <- if (is.na(pop_file)) pop_id else pop_file
  pop_id <- if (is.na(pop_id)) pop_file else pop_id
  sim_file <- if (is.na(sim_file)) sim_id else sim_file
  sim_id <- if (is.na(sim_id)) sim_file else sim_id

  # identify molecule column
  pop_mol <- ident_column(data, "pop.*(mol|molecule|molecules)$")
  sim_mol <- ident_column(data, "sim.*(mol|molecule|molecules)$")

  rep.fn(pop_id, "Population id")
  rep.fn(pop_file, "Population file")
  rep.fn(sim_id, "Simulation id")
  rep.fn(sim_file, "Simulation file")
  rep.fn(pop_mol, "Population molecules")
  rep.fn(sim_mol, "Simulation molecules")
  ############################################################################

  # Dose and Obs_Ids
  ############################################################################
  obs_ids <- ident_column(data, "obs.*(id|ids|identifier)$")
  rep.fn(obs_ids, "Observed identifiers")

  dose <- ident_column(data, "dose")
  rep.fn(dose, "Dose")
  dose.unit <- "mg"
  if (!is.na(dose)) {
    message("Extracting dosing unit ... ")
    tmp <- extract_unit(colnames(data)[dose])
    if (!has_units(tmp))
      stop("Could not extract dosing unit", call. = FALSE)

    dose.unit <- units::deparse_unit(tmp)
    if (!is_dose_unit(dose.unit))
      stop(paste("<", dose.unit ,"> is not a valid dosing unit"), call. = FALSE)
    else
      message(paste("<", dose.unit ,"> identified as dosing unit"))
  }

  # Plotting stuff
  ############################################################################
  header_col <- ident_column(data, "(.*)(main|header|headline|caption)")
  rep.fn(header_col, "Plot Header")

  x_min <- ident_column(data, "(x|time).*(min|from)")
  x_max <- ident_column(data, "(x|time).*(max|to)")
  x_unit <- ident_column(data, "(x|time).*(unit)")
  rep.fn(x_min, "Plot Minimum X")
  rep.fn(x_max, "Plot Maximum X")
  rep.fn(x_unit, "Plot X-Range Unit")

  y_min <- ident_column(data, "(y|value).*(min|from)")
  y_max <- ident_column(data, "(y|value).*(max|to)")
  y_unit <- ident_column(data, "(y|value).*(unit)")
  rep.fn(y_min, "Plot Minimum Y")
  rep.fn(y_max, "Plot Maximum Y")
  rep.fn(y_unit, "Plot Y-Range Unit")


  ############################################################################
  # Error checks
  n.na <- length(which(is.na(c(x_min, x_max, x_unit))))
  if (n.na > 0 && n.na < 3)
    stop("Plot Min/Max/Unit for X must be all defined or all undefined", call. = FALSE)

  n.na <- length(which(is.na(c(y_min, y_max, y_unit))))
  if (n.na > 0 && n.na < 3)
    stop("Plot Min/Max/Unit for Y must be all defined or all undefined", call. = FALSE)

  ############################################################################

  # Groups
  ############################################################################
  group_columns <- ident_column(data, "group", multi.matches = T)
  if (any_na(group_columns)) {
    message("No group columns identified")
    has_groups <- FALSE
  }
  else {
    message(paste("Identified group columns:", paste(group_columns, collapse = ", ")))
    has_groups <- TRUE
  }
  ############################################################################

  master <- master(id,
                   molecules = molecules,
                   reference = ref,
                   dose.unit = dose.unit,
                   groups = if (any_na(group_columns)) 0 else length(group_columns))

  for (row in 1:nrow(data)) {

    slice <- data %>% dplyr::slice(row)

    tryCatch({
      if (is.na(pop_mol)) {
        pop_molecules <- NA
      } else {
        entry <- trimws(slice[[pop_mol]])
        if (length(entry) == 0 || is.na(entry))
          pop_molecules <- NA
        else {
          str <- trimws(unlist(strsplit(entry, entry.sep, fixed = TRUE)))
          pop_molecules <- molecule_from_ids(molecules, str)
          if (length(pop_molecules) == 0)
            stop(paste("Could not find defined molecule with id <", str ,">\n"), call. = FALSE)
        }
      }

      if (is.na(sim_mol)) {
        sim_molecules <- NA
      } else {
        entry <- trimws(slice[[sim_mol]])
        if (length(entry) == 0 || is.na(entry))
          sim_molecules <- NA
        else {
          str <- trimws(unlist(strsplit(entry, entry.sep, fixed = TRUE)))
          sim_molecules <- molecule_from_ids(molecules, str)
          if (length(sim_molecules) == 0)
            stop(paste("Could not find defined molecule with id <", str ,">\n"), call. = FALSE)
        }
      }

      dosing <- if (is.na(dose)) NA else units::as_units(slice[[dose]], dose.unit)
      obs <- if (is.na(obs_ids)) NA else trimws(unlist(strsplit(slice[[obs_ids]],
                                                               entry.sep, fixed = TRUE)))


      plot.x.range <- NA
      if (!is.na(x_min))
        plot.x.range <- to_range(slice[[x_min]],
                                 slice[[x_max]],
                                 trimws(slice[[x_unit]]))

      plot.y.range <- NA
      if (!is.na(y_min))
        plot.y.range <- to_range(slice[[y_min]],
                                 slice[[y_max]],
                                 trimws(slice[[y_unit]]))

      groups <- NA
      if (has_groups) {
        groups <- slice[group_columns] %>% unlist(., use.names = FALSE) %>% trimws()
      }

      master <- add_master_entry(master.obj = master,
                                 pop.file = if (is.na(pop_file)) NA else trimws(slice[[pop_file]]),
                                 pop.id = if (is.na(pop_id)) NA else trimws(slice[[pop_id]]),
                                 pop.molecules = pop_molecules,
                                 sim.file = if (is.na(sim_file)) NA else trimws(slice[[sim_file]]),
                                 sim.id = if (is.na(sim_id)) NA else trimws(slice[[sim_id]]),
                                 sim.molecules = sim_molecules,
                                 dose = dosing,
                                 obs.ids = obs,
                                 groups = groups,
                                 plot.header = if (is.na(header_col)) NA else trimws(slice[[header_col]]),
                                 plot.x.range = plot.x.range,
                                 plot.y.range = plot.y.range

      )
    },
    error = function(e) {
      stop(paste("Error in row ", row, ":", e))
    })
  }

  return(master)
}



