#' Create a molecule object
#'
#' Constructs a `molecule` object containing metadata about a chemical compound,
#' optionally fetching the molecular weight and PubChem ID if not provided.
#'
#' @param name Character. Name of the molecule (mandatory, non-empty).
#' @param display.name Character. Display name of the molecule (defaults to `name` if `NA` or empty).
#' @param id Character. Identifier of the molecule (defaults to `name` if `NA` or empty).
#' @param file.name.match Character. Pattern used to match files (defaults to `name` if `NA` or empty).
#' This is used to identify the molecule file columns.
#' @param add.file.matcher Character vector. Additional file matching patterns.
#' This is used to identify the molecule file columns.
#' @param pubchem.id Numeric or character. PubChem compound ID (optional).
#' @param MW Numeric. Molecular weight (optional; fetched from PubChem if `NA`).
#' @param is.fraction Logical. Whether the molecule is expressed as a fraction (default `FALSE`).
#' E.g. for fraction extrected in urin.
#' @param fixed.unit `units` object. Fixed unit of measurement (optional; must be of class `units` if provided).
#' @param ylab Character. Label for y-axis (default `"Plasma Concentration"`).
#' @param color Character. Colour for plotting (default `"black"`).
#' @param pch Integer. Point character for plotting (default `19`).
#' @param lty Integer. Line type for plotting (default `1`).
#' @param in.legend Logical. Whether to include this molecule in plot legends (default `TRUE`).
#'
#' @return An object of class `molecule`, which is a list with elements:
#' \itemize{
#' \item name
#' \item display.name
#' \item id
#' \item file.name.match
#' \item add.file.matcher
#' \item pubchem.id
#' \item MW
#' \item fixed.unit
#' \item is.fraction
#' \item color
#' \item ylab
#' \item pch
#' \item lty
#' \item in.legend
#' }
#'
#' @details
#' If `MW` is not provided, the function queries the PubChem REST API using either
#' `pubchem.id` or `name` to retrieve the molecular weight and PubChem ID.
#' If multiple PubChem entries are found, the function stops with an error.
#'
#' @examples
#' # Minimal example
#' \dontrun{
#' mol <- molecule("Caffeine")
#'
#' # With fixed MW and display name
#' mol <- molecule("Caffeine", display.name = "Caffeine (custom)", MW = 194.19)
#' }
#' @export
molecule <- function(name,
                     display.name = NA,
                     id = NA,
                     file.name.match = NA,
                     add.file.matcher = c(),
                     pubchem.id = NA,
                     MW = NA,
                     is.fraction = FALSE,
                     fixed.unit = NA,
                     ylab = "Plasma Concentration",
                     color = "black",
                     pch = 19,
                     lty = 1,
                     in.legend = TRUE) {
  name <- paste(trimws(name))

  add.file.matcher <- c(add.file.matcher)

  if (nchar(trimws(name)) == 0) {
    stop("Molecule needs a non-empty name")
  }

  if (is.na(display.name) || nchar(trimws(display.name)) == 0) {
    display.name <- name
  } else {
    display.name <- paste(trimws(display.name))
  }

  if (is.na(file.name.match) || nchar(trimws(file.name.match)) == 0) {
    file.name.match <- name
  } else {
    file.name.match <- paste(trimws(file.name.match))
  }

  if (is.na(id) || nchar(trimws(id)) == 0) {
    id <- name
  } else {
    id <- paste(trimws(id))
  }

  if (is.na(MW)) {
    message("Fetching data from Pubchem")
    json.data <- NA
    req <- "/property/MolecularWeight/JSON"
    req.stem <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
    if (!is.na(pubchem.id)) {
      json.data <- httr::GET(paste(req.stem, "cid/", pubchem.id, req, sep = ""))
    } else {
      json.data <- httr::GET(paste(req.stem, "name/", name, req, sep = ""))
    }

    data <- httr::content(json.data, "text", encoding = "UTF-8")
    data <- as.data.frame(jsonlite::fromJSON(data, flatten = TRUE))

    pubchem.id <- data$PropertyTable.Properties.CID
    MW <- data$PropertyTable.Properties.MolecularWeight

    if (length(pubchem.id) > 1) {
      stop(paste(
        "Fetching data from Pubchem for molecule <", name, ">: Ambiguous entries found with ids: ",
        paste(pubchem.id, collapse = ", ")
      ), call. = F)
    }

    if (is.null(MW)) {
      stop(paste("Fetching data from Pubchem for molecule <", name, "> failed"))
    }
  }

  if (is.na(MW) || !is.numeric(MW)) {
    stop("MW must be numeric and not NA")
  }


  if (!is.logical(is.fraction)) {
    stop("is.fraction must be logical")
  }

  if (!is.logical(in.legend)) {
    stop("in.legend must be logical")
  }

  if (!is.na(fixed.unit) && !.has.units(fixed.unit)) {
    stop("fixed.unit must be of class units. Did you forget as_units(...)?")
  }

  result <- list(
    name = name,
    display.name = display.name,
    id = id,
    file.name.match = file.name.match,
    add.file.matcher = add.file.matcher,
    pubchem.id = pubchem.id,
    color = color,
    ylab = ylab,
    pch = pch,
    lty = lty,
    MW = MW,
    fixed.unit = fixed.unit,
    is.fraction = is.fraction,
    in.legend = in.legend
  )

  class(result) <- "molecule"
  return(result)
}

#' Check if an Object is a Molecule
#'
#' This function checks whether the given object inherits from the class "molecule".
#'
#' @param x An object to be checked.
#' @return A logical value: `TRUE` if the object inherits from the class "molecule", otherwise `FALSE`.
#' @export
is.molecule <- function(x) inherits(x, "molecule")


#' Clone a Molecule Definition
#'
#' Creates a copy of an existing molecule definition allowing overriding of
#' selected fields.
#'
#' @param molecule The molecule to copy.
#' @param name Optional new name.
#' @param display.name Optional new display name.
#' @param id Optional new identifier.
#' @param file.name.match Optional new file matching pattern.
#' @param add.file.matcher Additional file matchers to append.
#' @param pubchem.id Optional PubChem ID.
#' @param MW Optional molecular weight.
#' @param is.fraction Whether the molecule represents a fraction.
#' @param fixed.unit Fixed unit for this molecule.
#' @param ylab Y-axis label used in plots.
#' @param color Plotting color.
#' @param pch Point type for plotting.
#' @param lty Line type for plotting.
#' @param in.legend Should the molecule appear in legends.
#'
#' @return A new molecule object.
#' @export
molecule_clone <- function(molecule,
                           name = NULL,
                           display.name = NULL,
                           id = NULL,
                           file.name.match = NULL,
                           add.file.matcher = NULL,
                           pubchem.id = NULL,
                           MW = NULL,
                           is.fraction = NULL,
                           fixed.unit = NULL,
                           ylab = NULL,
                           color = NULL,
                           pch = NULL,
                           lty = NULL,
                           in.legend = NULL) {
  mol <- molecule

  if (!is.null(name)) {
    mol$name <- paste(trimws(name))
  }

  if (!is.null(id)) {
    mol$id <- paste(trimws(id))
  }

  if (!is.null(name)) {
    mol$name <- name
  }

  if (!is.null(file.name.match)) {
    mol$file.name.match <- paste(trimws(file.name.match))
  }

  if (!is.null(add.file.matcher)) {
    mol$add.file.matcher <- c(add.file.matcher)
  }

  if (!is.null(pubchem.id)) {
    mol$pubchem.id <- pubchem.id
  }

  if (!is.null(MW)) {
    if (is.na(MW) || !is.numeric(MW)) {
      stop("MW must be numeric and not NA")
    }
    mol$MW <- MW
  }

  if (!is.null(is.fraction)) {
    if (!is.logical(is.fraction)) {
      stop("is.fraction must be logical")
    }

    mol$is.fraction <- is.fraction
  }

  if (!is.null(fixed.unit)) {
    mol$fixed.unit <- fixed.unit
  }

  if (!is.null(ylab)) {
    mol$ylab <- ylab
  }

  if (!is.null(color)) {
    mol$color <- color
  }

  if (!is.null(pch)) {
    mol$pch <- pch
  }

  if (!is.null(lty)) {
    mol$lty <- lty
  }

  if (!is.null(in.legend)) {
    if (!is.logical(in.legend)) {
      stop("in.legend must be logical")
    }

    mol$in.legend <- in.legend
  }

  return(mol)
}


#' Combine Matched Profile Objects
#'
#' Utility to merge several `MatchedProfiles` objects by group. The function can
#' optionally rename IDs and choose a reference profile.
#'
#' @param profiles List of `MatchedProfiles` objects to combine.
#' @param by_group Group index used to determine which profiles belong together.
#' @param rename_id_by_group Optional group index whose value becomes the new ID.
#' @param ref_group If provided, selects the profile with `ref_id` from this
#'   group as template.
#' @param ref_id Name of the reference group.
#' @param silent Logical controlling progress messages.
#'
#' @return A list of combined `MatchedProfiles` objects.
#' @export
combine_profiles <- function(profiles,
                             by_group = 1,
                             rename_id_by_group = NULL,
                             ref_group = NULL,
                             ref_id = "control",
                             silent = FALSE) {
  # prepare
  msg_fn <- if (silent) {
    function(x, nl = TRUE) {}
  } else {
    function(x, nl = TRUE) {
      message(x, appendLF = nl)
    }
  }

  # error checks
  if (!is.numeric(by_group) || by_group <= 0) {
    stop("by_group must be positiv numeric (e.g. 1 for group 1", call. = FALSE)
  }

  if (!is.null(ref_group) && (!is.numeric(by_group) || by_group <= 0)) {
    stop("ref_group must be positiv numeric (e.g. 1 for group 1", call. = FALSE)
  }

  # gather groups
  group_ids <- c()
  for (pro in profiles) {
    if (length(pro$groups) < by_group) {
      stop("by_group idx is > lenght of profile groups", call. = FALSE)
    }

    group_ids <- c(group_ids, pro$groups[by_group])
  }

  new_profiles <- list()
  groups <- unique(group_ids)
  msg_fn(paste("Found <", length(groups), "> groups"))

  # iterate over all combine-groups
  for (g in groups) {
    msg_fn(paste("* Combine profiles of group <", g, ">"))

    pro_list <- list()
    idxs <- which(group_ids == g)
    for (idx in idxs) {
      pro_list <- append(pro_list, list(profiles[[idx]]))
    }

    # check for ref
    if (is.null(ref_group)) {
      combined_pro <- pro_list[[1]]
      pro_list[[1]] <- NULL
    } else {
      combined_pro <- NULL
      for (i in seq_along(pro_list)) {
        if (length(pro_list[[i]]$groups) < ref_group) {
          stop("ref_group idx is > lenght of profile groups", call. = FALSE)
        }

        if (pro_list[[i]]$groups[ref_group] == ref_id) {
          combined_pro <- pro_list[[i]]
          pro_list[[i]] <- NULL
          break
        }
      }

      if (is.null(combined_pro)) {
        stop(paste("Could not find ref_id <", ref_id, "> in profiles of group", call. = FALSE))
      }
    }

    # we gather in the first profile (profiles + observed ids)
    x.offset <- list(
      profiles = rep(combined_pro$plot.infos$x.offset, length(combined_pro$profiles)),
      obs = c(rep(combined_pro$plot.infos$x.offset, length(combined_pro$obs.ids)))
    )
    for (pro in pro_list) {
      combined_pro$profiles <- append(combined_pro$profiles, pro$profiles)

      x.offset$profiles <- c(x.offset$profiles, rep(
        pro$plot.infos$x.offset,
        length(pro$profiles)
      ))

      if (length(pro$obs.ids) > 0) {
        combined_pro$obs.ids <- append(combined_pro$obs.ids, pro$obs.ids)

        x.offset$obs <- c(x.offset$obs, rep(
          pro$plot.infos$x.offset,
          length(pro$obs.ids)
        ))
      }
    }
    combined_pro$plot.infos$x.offset <- x.offset
    combined_pro$plot.infos$x.offset$obs[is.na(combined_pro$plot.infos$x.offset$obs)] <- 0
    combined_pro$plot.infos$x.offset$profiles[is.na(combined_pro$plot.infos$x.offset$profiles)] <- 0

    # rename id
    if (!is.null(rename_id_by_group)) {
      if (!is.numeric(rename_id_by_group) || rename_id_by_group <= 0) {
        stop("rename_id_by_group must be positiv numeric (e.g. 1 for group 1", call. = FALSE)
      }


      if (length(combined_pro$groups) < rename_id_by_group) {
        stop("rename_id_by_group idx is > lenght of profile group", call. = FALSE)
      }

      combined_pro$id <- combined_pro$groups[rename_id_by_group]
    }

    # Check if we have duplicated OBS-IDs
    if (any(duplicated(combined_pro$obs.ids))) {
      stop("Combined profiles with the same observed data IDs", call. = FALSE)
    }


    new_profiles <- append(new_profiles, list(combined_pro))
  }

  return(new_profiles)
}
