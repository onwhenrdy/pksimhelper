


#' Title
#'
#' @param name
#' @param display.name
#' @param id
#' @param file.name.match
#' @param add.file.matcher
#' @param pubchem.id
#' @param MW
#' @param is.fraction
#' @param fixed.unit
#' @param ylab
#' @param color
#' @param pch
#' @param lty
#' @param in.legend
#'
#' @return
#' @export
#'
#' @examples
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

  if (nchar(trimws(name)) == 0)
    stop("Molecule needs a non-empty name")

  if (is.na(display.name) || nchar(trimws(display.name)) == 0)
    display.name <- name
  else
    display.name <- paste(trimws(display.name))

  if (is.na(file.name.match) || nchar(trimws(file.name.match)) == 0)
    file.name.match <- name
  else
    file.name.match <- paste(trimws(file.name.match))

  if (is.na(id) || nchar(trimws(id)) == 0)
    id <- name
  else
    id <- paste(trimws(id))

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
      stop(paste("Fetching data from Pubchem for molecule <", name, ">: Ambiguous entries found with ids: ",
                 paste(pubchem.id, collapse = ", ")), call. = F)
    }

    if (is.null(MW))
      stop(paste("Fetching data from Pubchem for molecule <", name, "> failed"))
  }

  if (is.na(MW) || !is.numeric(MW))
    stop("MW must be numeric and not NA")


  if(!is.logical(is.fraction)) {
    stop("is.fraction must be logical")
  }

  if(!is.logical(in.legend)) {
    stop("in.legend must be logical")
  }

  if(!is.na(fixed.unit) && !.has.units(fixed.unit)) {
    stop("fixed.unit must be of class units. Did you forget as_units(...)?")
  }

  result <- list(name = name,
                 display.name = display.name,
                 id = id,
                 file.name.match = file.name.match,
                 add.file.matcher =  add.file.matcher,
                 pubchem.id = pubchem.id,
                 color = color,
                 ylab = ylab,
                 pch = pch,
                 lty = lty,
                 MW = MW,
                 fixed.unit = fixed.unit,
                 is.fraction = is.fraction,
                 in.legend = in.legend)

  class(result) <- "molecule"
  return(result)
}

# molecule helper
is.molecule <- function(x) inherits(x, "molecule")


# clones a molecule
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

    if (!is.null(name))
      mol$name <- paste(trimws(name))

    if (!is.null(id))
      mol$id <- paste(trimws(id))

    if (!is.null(name))
      mol$name <- name

    if (!is.null(file.name.match))
      mol$file.name.match <- paste(trimws(file.name.match))

    if (!is.null(add.file.matcher))
      mol$add.file.matcher <- c(add.file.matcher)

    if (!is.null(pubchem.id))
      mol$pubchem.id <- pubchem.id

    if (!is.null(MW)) {
      if (is.na(MW) || !is.numeric(MW))
        stop("MW must be numeric and not NA")
      mol$MW <- MW
    }

    if (!is.null(is.fraction)) {
      if(!is.logical(is.fraction))
        stop("is.fraction must be logical")

      mol$is.fraction <- is.fraction
    }

    if (!is.null(fixed.unit))
      mol$fixed.unit <- fixed.unit

    if (!is.null(ylab))
      mol$ylab <- ylab

    if (!is.null(color))
      mol$color <- color

    if (!is.null(pch))
      mol$pch <- pch

    if (!is.null(lty))
      mol$lty <- lty

    if (!is.null(in.legend)) {
      if(!is.logical(in.legend))
        stop("in.legend must be logical")

      mol$in.legend <- in.legend
    }

    return(mol)
}


combine_profiles <- function(profiles,
                             by_group = 1,
                             rename_id_by_group = NULL,
                             ref_group = NULL,
                             ref_id = "control",
                             silent = FALSE) {

  # prepare
  msg_fn <- if (silent) function(x, nl = TRUE) {} else function(x, nl = TRUE) {message(x, appendLF = nl)}

  # error checks
  if (!is.numeric(by_group) || by_group <= 0)
    stop("by_group must be positiv numeric (e.g. 1 for group 1", call.= FALSE)

  if (!is.null(ref_group) && (!is.numeric(by_group) || by_group <= 0))
    stop("ref_group must be positiv numeric (e.g. 1 for group 1", call.= FALSE)

  # gather groups
  group_ids <- c()
  for (pro in profiles) {

    if (length(pro$groups) < by_group)
      stop("by_group idx is > lenght of profile groups", call. = FALSE)

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
    } else  {

      combined_pro <- NULL
      for (i in seq_along(pro_list)) {
        if (length(pro_list[[i]]$groups) < ref_group)
          stop("ref_group idx is > lenght of profile groups", call. = FALSE)

        if (pro_list[[i]]$groups[ref_group] == ref_id) {
          combined_pro <- pro_list[[i]]
          pro_list[[i]] <- NULL
          break
        }
      }

      if (is.null(combined_pro))
        stop(paste("Could not find ref_id <", ref_id, "> in profiles of group", call. = FALSE))
    }

    # we gather in the first profile (profiles + observed ids)
    x.offset <- list(profiles = rep(combined_pro$plot.infos$x.offset, length(combined_pro$profiles)),
                     obs = c(rep(combined_pro$plot.infos$x.offset, length(combined_pro$obs.ids))))
    for (pro in pro_list) {

      combined_pro$profiles <- append(combined_pro$profiles, pro$profiles)

      x.offset$profiles <- c(x.offset$profiles, rep(pro$plot.infos$x.offset,
                                                    length(pro$profiles)))

      if(!is.na(pro$obs.ids)) {
        combined_pro$obs.ids <- append(combined_pro$obs.ids, pro$obs.ids)

        x.offset$obs <- c(x.offset$obs, rep(pro$plot.infos$x.offset,
                                            length(pro$obs.ids)))
      }
    }
    combined_pro$plot.infos$x.offset <- x.offset
    combined_pro$plot.infos$x.offset$obs[is.na(combined_pro$plot.infos$x.offset$obs)] <- 0
    combined_pro$plot.infos$x.offset$profiles[is.na(combined_pro$plot.infos$x.offset$profiles)] <- 0

    # rename id
    if(!is.null(rename_id_by_group)) {

      if (!is.numeric(rename_id_by_group) || rename_id_by_group <= 0)
        stop("rename_id_by_group must be positiv numeric (e.g. 1 for group 1", call.= FALSE)


      if (length(combined_pro$groups) < rename_id_by_group)
        stop("rename_id_by_group idx is > lenght of profile group", call. = FALSE)

      combined_pro$id <- combined_pro$groups[rename_id_by_group]
    }

    # Check if we have duplicated OBS-IDs
    if (any(duplicated(combined_pro$obs.ids)))
      stop("Combined profiles with the same observed data IDs", call. = FALSE)


    new_profiles <- append(new_profiles, list(combined_pro))
  }

  return(new_profiles)
}




