
#################################################################################
# Generic Helper
#################################################################################

#' Fetches Pubchem Information from ONE Molecule Name or CID
#'
#' @param name_or_cid Name of CID of a molecule.
#' @param features The desired features (empty will only fetch CID)
#'
#' @return A `tibble` of the desired features including the CID as the first entry.
#' @noRd
#'
.fetch_properties_from_pubchem <- function(name_or_cid,
                                features = c("MolecularWeight")) {

  if (length(name_or_cid) != 1)
    stop("fetch_from_pubchem will only work for one query (not vectorized)", call. = F)

  if (is.na(name_or_cid[1]))
    return(NA)

  # remove NA from features
  features <- features[!is.na(features)]

  if (is.na(features[1]) || is.null(features[1]))
    feature_str <- ""
  else
    feature_str <-  paste(features, collapse = ",")


  req <- if (nchar(feature_str) == 0) "/cids/JSON" else paste0("/property/", feature_str, "/JSON")
  req.stem <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"

  json.data <- NA
  if (is.numeric(name_or_cid)) {
    json.data <- httr::GET(paste0(req.stem, "cid/", name_or_cid, req))
  } else {
    json.data <- httr::GET(paste0(req.stem, "name/", utils::URLencode(name_or_cid), req))
  }

  data <- httr::content(json.data, "text", encoding = "UTF-8")
  data <- as.data.frame(jsonlite::fromJSON(data, flatten = TRUE))
  # error handling
  if (!is.null(data$Fault.Message))
    stop(paste("Fetching data from Pubchem for <",
               name_or_cid, "> failed:", data$Fault.Message), call. = F)

  # if no CID is found the API will just return the CID entry
  if (ncol(data) != length(features) + 1)
    stop(paste("Fetching data from Pubchem for <",
               name_or_cid, "> failed: Could not find entry."), call. = F)


  colnames(data) <- sub("PropertyTable.Properties.", "", colnames(data))
  return(tibble::as_tibble(data))
}

#' Fetches Compound properties from PubChem
#'
#' Fetches properties like e.g. MW or XlogP3 from Pubchem for compounds or CIDs.
#'
#' @param names_or_cids Compound name or CID or vector/list-like compound names and/or CIDs
#' @param features The requested features from PubChem.
#'   For a list of possbile features see \url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$_Toc494865567}.
#'   If no features are provided `NA`, `NULL` or `c` only the CIDs are returned.
#'
#' @return `tibble` with first column CID and requested feature columns.
#' @export
#'
#' @examples
#' nic_data <- fetch_properties_from_pubchem("Nicotine")
#' mol_data <- fetch_properties_from_pubchem(list("Nicotine", 12345),
#'                                           features = c("MolecularWeight", "XLogP"))
#'
fetch_properties_from_pubchem <- function(names_or_cids, features = c("MolecularWeight")) {
  message("Fetching data from Pubchem")
  data <- lapply(names_or_cids, .fetch_properties_from_pubchem, features = features)
  do.call(rbind, data)
}


#################################################################################
# Molecule
#################################################################################

#' Creates a Molecule
#'
#' A molecule is a description of a measured entity from a PKSim/Mobi container.
#'   It can be a PK descriptor (e.g. concentrations) or PD descriptor (e.g. heart rate).
#'
#' If MW is not set, the name or (if set) the pubchem identifier will be
#'   used to fetch the molecular weight from the Pubchem database. If this fails the molecule is
#'   not created. If the fetching is done via name (no pubchem identifier is provided) the pubchem
#'   idenfier is fetched as well.
#'
#' For PD Molecules a dummy value for MW should be set (e.g. 0) to avoid the automatic fetching
#'   from Pubchem. When it comes to matchers, special charaters (e.g. "|") must be regular
#'   expression safe (escaped). We povide special matcher functions that allows the easy creation
#'   often used matching expressions.
#'
#' @param name The name of the molecule (must be a single non-empty string).
#' @param display.name The name that should be displayed in plots (e.g. in the legend). Must be
#'   a single non-empty string.
#' @param id A unique id for the Molecule. Must be a single value that will be converted to string.
#' @param container A single non-empty string or NA that defines the model container the molecule
#'   is simulated.
#' @param column.matcher Strings that match the Molecule for the parser functions (AND relationship).
#'   Must be a vector/single value of non-empty strings.
#' @param pubchem.id The pubchem identifier. Must be a single numeric value or `NA`.
#' @param MW The molecular weight of the Molecule. Must be a single numeric value or `NA`.
#' @param is.fraction Single logical value that indicates if the Molecule represents a fraction
#'   (e.g. fraction excreted to urin).
#' @param fixed.unit Indicates if the Molecule has a fixed unit (always the same unit).
#'   If not a single unit is used `NA` must be set.
#' @param ylab A label for the y-axis if Profiles for a Molecule are plotted. Must be a single
#'   non-empty string.
#' @param col A color that should be used for plotting. Must be a single value (`NA` is allowed).
#' @param pch A point type that should be used for plotting. Must be a single value (`NA` is allowed).
#' @param lty A line type that should be used for plotting. Must be a single value (`NA` is allowed).
#'
#' @export
#'
#' @family molecule functions
#'
#' @examples
#' nic <- molecule("Nicotine") # will fetch from PubChem
#' nic_2 <- molecule("Nic", pubchem.id = 89594) # will fetch Nicotine from Pubchem from id
#' test <- molecule("PD Marker", MW = 0) # will not fetch from Pubchem
#'
molecule <- function(name,
                     display.name = name,
                     id = name,
                     container = NA,
                     column.matcher = c(name),
                     pubchem.id = NA,
                     MW = NA,
                     is.fraction = F,
                     fixed.unit = NA,
                     ylab = "Plasma Concentration",
                     col = "black",
                     pch = 19,
                     lty = 1) {

  name <- trimws(name)
  if (!is_single_string(name, can.be.empty = F))
    stop("name must be a non-empty single string", call. = F)

  display.name <- trimws(display.name)
  if (!is_single_string(display.name, can.be.empty = F))
    stop("display.name must be a non-empty single string", call. = F)

  id <- trimws(id)
  if (!is_single_string(id, can.be.empty = F))
    stop("id must be a non-empty", call. = F)

  container <- trimws(container)
  if (!is_single_string(id, can.be.empty = F) && !any_na(container))
    stop("container must be a non-empty string or NA", call. = F)

  if (length(column.matcher) == 0)
    stop("column.matcher must be non-empty", call. = F)
  column.matcher <- trimws(column.matcher)
  column.matcher <- as.vector(column.matcher)
  if (!all(sapply(column.matcher, is_single_string, can.be.empty = F)))
    stop("column.matcher cannot have non-empty entries", call. = F)

  if (!is_single_numeric(pubchem.id) && !any_na(pubchem.id))
    stop("pubchem.id must be a single numeric value or NA", call. = F)

  if (length(fixed.unit) != 1 || (!has_units(fixed.unit) && !any_na(fixed.unit)))
    stop("fixed.unit must be a single unit or NA", call. = F)

  if (!is_single_numeric(MW) && !any_na(MW))
    stop("MW must be a single numeric value or NA", call. = F)

  if (!is_single_logical(is.fraction))
    stop("is.fraction must be a single logical value", call. = F)

  ylab <- trimws(ylab)
  if (!is_single_string(ylab, can.be.empty = F))
    stop("ylab must be a non-empty single string", call. = F)

  if (!is_single_string(ylab, can.be.empty = F))
    stop("ylab must be a non-empty single string", call. = F)

  # NAs are allowed
  if (length(pch) != 1 || !is.vector(pch))
    stop("pch must be a single value", call. = F)

  if (length(lty) != 1 || !is.vector(pch))
    stop("lty must be a single value", call. = F)

  if (length(col) != 1 || !is.vector(col))
    stop("col must be a single value", call. = F)

  if (is.na(MW)) {
    query <- if (!any_na(pubchem.id)) pubchem.id else name
    dat <- fetch_properties_from_pubchem(query)
    pubchem.id <- dat$CID
    MW <- dat$MolecularWeight
  }

  result <- list(name = name,
                 display.name = display.name,
                 id = id,
                 container = container,
                 column.matcher = column.matcher,
                 pubchem.id = pubchem.id,
                 col = col,
                 ylab = ylab,
                 pch = pch,
                 lty = lty,
                 MW = MW,
                 fixed.unit = fixed.unit,
                 is.fraction = is.fraction)

  class(result) <- append("molecule", class(result))
  return(result)
}


#' Prints Information about a Molecule
#'
#' @param x The molecule object.
#' @param details True if all details should be printed, else False.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisible copy of the object \code{x}.
#'
#' @family molecule functions
#'
#' @export
#'
#' @examples
#' a <- molecule("Nicotine")
#' a # prints a without details
#' print(a, details = TRUE) # prints a with details
#'
#' b <- molecule("THC")
#' mol_list <- list(a, b)
#' print(mol_list, details = TRUE) # prints molecule list with details
#'
print.molecule <- function(x, details = F, ...) {

  if (!is_molecule(x))
    stop("Input must be of class molecule", call. = F)

  if (!is_single_logical(details))
    stop("details must be logical", call. = F)

  cat(paste(x$name, "(Molecule)\n"))
  cat(paste("Id:\t\t", x$id, "\n"))
  cat(paste("MW:\t\t", x$MW, "\n"))
  cat(paste("Container:\t", x$container))

  if (details) {
    cat(paste("\nDisplay:\t", x$display.name, "\n"))
    cat(paste("Matches:\t", paste(x$column.matcher, collapse = ", "), "\n"))
    cat(paste("Pubchem ID:\t", x$pubchem.id, "\n"))
    cat(paste("Color:\t\t", x$col, "\n"))
    cat(paste("Y label:\t", x$ylab, "\n"))
    cat(paste("Point type:\t", x$pch, "\n"))
    cat(paste("Line type:\t", x$lty, "\n"))
    cat(paste("Fixed unit:\t", x$fixed.unit, "\n"))
    cat(paste("Fraction:\t", x$is.fraction))
  }

  invisible(x)
}


#' Checks if a Data Structure is a Molecule
#'
#' @param ds The data structure that should be tested.
#'
#' @return `TRUE` if the data structure is a Molecule.
#' @export
#'
#' @family molecule functions
#'
#' @examples
#' is_molecule("Hallo") # FALSE
#' is_molecule(molecule("Test", MW = 100)) # TRUE
#'
is_molecule <- function(ds) {
  inherits(ds, "molecule")
}


#' Creates a List of Molecules
#'
#' @param ... Must be a non-empty list of molecules.
#'
#' @return A list of molecules that are named by id
#' @export
#'
#' @family molecule functions
#'
#' @examples
#'
#' mol_1 <- molecule("Mol 1", MW = 100)
#' mol_2 <- molecule("Mol 2", MW = 200)
#' mol_list <- molecule_list(mol_1, mol_2)
#' is_molecule_list(mol_list) # is TRUE
#'
molecule_list <- function(...) {

  mols <- list(...)
  if (!is_molecule_list(mols))
    stop("arguments are not of class molecule", call. = F)

  names(mols) <- sapply(mols, function(x) x$id)
  return(mols)
}

#' Returns a Vector of Molecule Ids from a Molecule List
#'
#' @param mol.list A non-empty molecule list
#'
#' @return A vector of molecule Ids.
#' @export
#'
#' @family molecule functions
#'
#' @examples
#'
#' mol_1 <- molecule("Mol 1", MW = 100)
#' mol_2 <- molecule("Mol 2", MW = 200)
#' mol_list <- molecule_list(mol_1, mol_2)
#' ids <- molecule_ids(mol_list)
#'
molecule_ids <- function(mol.list) {

  if (!is_molecule_list(mol.list))
    stop("mol_list must be a molecule list", call. = F)

  unname(sapply(mol.list, function(x) x$id))
}


#' Returns Molecules with certain ids from a Molecule List
#'
#' @param mol.list A non-empty molecule list.
#' @param query.ids One or more ids that should be queried
#' @param partial.matches If partial matches are allowed (will return `NULL` entries for non-matches).
#' @param rm.na If non-matches should be removed from the result (only works for
#'   `partial.matches = TRUE`).
#'
#' @return If `query_ids` is one string a molecule is returned. If more than one molecule is queried
#'   a molecule list is returned. If one or more query strings could not be found `NULL` is returned
#'   (`partial.matches = FALSE`).
#' @export
#'
#' @family molecule functions
#'
#' @examples
#'
#' mol_1 <- molecule("Mol 1", MW = 100)
#' mol_2 <- molecule("Mol 2", MW = 200)
#' mol_list <- molecule_list(mol_1, mol_2)
#'
#' molecule_from_ids(mol_list, c("Mol 1")) # only one molecule is returned
#' molecule_from_ids(mol_list, c("Mol 1", "Mol 2")) # two molecules are returned
#' molecule_from_ids(mol_list, c("Mol 3", "Mol 2")) # NULL is returned
#'
#' # list(NULL, Mol_2) is returned
#' molecule_from_ids(mol_list, c("Mol 3", "Mol 2"), partial.matches = TRUE)
#'
#' # Mol_2 is returned
#' molecule_from_ids(mol_list, c("Mol 3", "Mol 2"), partial.matches = TRUE, rm.na = TRUE)
#'
molecule_from_ids <- function(mol.list, query.ids, partial.matches = FALSE, rm.na = FALSE) {

  if (!is_string_list(query.ids, can.be.empty = T) || length(query.ids) == 0)
    stop("query_ids must be a single value or vector of non-empty string", call. = F)

  if (!is_single_logical(partial.matches))
    stop("partial.matches must be a single logical value", call. = F)

  if (!is_single_logical(rm.na))
    stop("rm.na must be a single logical value", call. = F)

  ids <- molecule_ids(mol.list)
  ids <- tolower(ids)
  query.ids <- tolower(query.ids)
  idx <- base::match(query.ids, ids)
  if (!partial.matches && any_na(idx))
    return(NULL)

  if (rm.na)
    idx <- idx[!is.na(idx)]

  if (length(idx) == 0)
    return(NULL)

  if (length(idx) == 1)
    return(mol.list[[idx]])

  return(mol.list[idx])
}


#' Checks if a Molecule List Contains all Query Ids
#'
#' @param mol.list A non-empty molecule list.
#' @param query.ids One or more ids that should be queried.
#'
#' @return `TRUE` if all query ids could be found, else `FALSE`.
#' @export
#'
#' @family molecule functions
#'
#' @examples
#'
#' mol_1 <- molecule("Mol 1", MW = 100)
#' mol_2 <- molecule("Mol 2", MW = 200)
#' mol_list <- molecule_list(mol_1, mol_2)
#' has_molecules(mol_list, c("Mol 1")) # TRUE
#' has_molecules(mol_list, c("Mol 1", "Mol 2")) # TRUE
#' has_molecules(mol_list, c("Mol 3", "Mol 2")) # FALSE
#'
has_molecules <- function(mol.list, query.ids) {

  length(molecule_from_ids(mol.list, query.ids)) != 0
}


#' Merges two Molecule Lists
#'
#' Trick: You can make a molecule list unique by merging the list with itself.
#'
#' @param list.a The first molecule list (must be non-empty).
#' @param list.b The second molecule list (must be non-empty).
#'
#' @return A merged molecule list with unique ids.
#' @export
#'
#' @family molecule functions
#'
#' @examples
#'
#' mol_1 <- molecule("Mol 1", MW = 100)
#' mol_2 <- molecule("Mol 2", MW = 200)
#' mol_list <- molecule_list(mol_1, mol_2)
#' merge_molecule_lists(mol_list, mol_list) # returns mol_list
#'
merge_molecule_lists <- function(list.a, list.b) {

  ids_a <- molecule_ids(list.a)
  ids_b <- molecule_ids(list.b)

  result <- c(list.a, list.b)
  duplicates <- base::duplicated(c(ids_a, ids_b))
  result <- result[!duplicates]

  return(result)
}

#' Checks if a Molecule List has Unique Ids
#'
#' @param mol.list The molecule list that should be tested.
#'
#' @return `TRUE`, if the molecule list is unique, else `FALSE`.
#' @export
#'
#' @family molecule functions
#'
is_unique_molecule_list <- function(mol.list) {

  ids <- molecule_ids(mol.list)
  return(length(ids) == length(unique(ids)))
}

#' Checks if a Data Structure is a list of Molecules
#'
#' @param list The data structure that should be tested.
#'
#' @return `TRUE` if the data structure is a Molecule list. For single Molecules
#'   this function returns `FALSE`.
#' @export
#'
#' @family molecule functions
#'
is_molecule_list <- function(list) {

  if (length(list) < 1)
    return(F)

  all(sapply(list, is_molecule))
}

#################################################################################
# Master data
#################################################################################

#' Creates a Master Object
#'
#' The reference should typically store the filename if the master object is parsed from a file.
#'
#' @param id The id of the master object (must be a non-empty string).
#' @param reference The reference of the master object (optional).
#' @param groups The number of Groups per master entry that should be created
#'   (must be numeric and at least be 0).
#' @param molecules A molecule list or single molecule that should at least contain all molecules
#'   that are added later on. Must be non-empty.
#' @param dose.unit A valid unit string for a given dose. Unit string must be lower case.

#' Adding entries to a master object should only be done via \code{add_master_entry}.
#'
#' @return A master object.
#' @export
#'
#' @family master functions
#'
#' @examples
#'
#' mol <- molecule("Mol 1", MW = 100)
#' m1 <- master("test", mol)
#' m2 <- master("test 2", mol, "my ref", groups = 3)
#' m3 <- master("test 2", mol, "my ref", "ng",  groups = 3)
#'
master <- function(id, molecules, reference = NA, dose.unit = "mg", groups = 0) {

  if (!is_single_string(id, can.be.empty = F)) {
    stop("id must be a single non-empty string", call. = F)
  }

  if (!is_dose_unit(dose.unit)) {
    stop("dose_unit must be a valid unit string that is convertible to mg", call. = F)
  }

  if (!is_single_numeric(groups) || groups < 0)
    stop("groups must be >= 0", call. = F)

  if (is_molecule(molecules))
    molecules <- molecule_list(molecules)

  if (!is_molecule_list(molecules))
    stop("molecules must be a non-empty molecule list", call. = F)

  if (!is_unique_molecule_list(molecules))
    stop("molecules must be unique - double entries found", call. = F)

  # basic names
  col.names <- c("Pop_Id", "Pop_File", "Pop_Mols", "Sim_Id", "Sim_File",
                 "Sim_Mols", "Dose", "Obs_Ids")

  group.names <- c()
  # groups
  if (groups > 0) {
    group.names <- paste0("Group_", seq(1, groups))
    col.names <- c(col.names, group.names)
  }

  # plot names
  col.names <- c(col.names, "Plot_Header", "Plot_X_Min", "Plot_X_Max", "Plot_X_Unit",
                 "Plot_Y_Min", "Plot_Y_Max", "Plot_Y_Unit")

  mat <- matrix(ncol = length(col.names), nrow = 0)
  colnames(mat) <- col.names
  df <- tibble::as_tibble(mat)

  result <- list(id = id,
                 dose.unit = dose.unit,
                 reference = reference,
                 group.names = group.names,
                 molecules = molecules,
                 data = df)

  class(result) <- append("master", class(result))
  return(result)
}

#' Generic Function to Return Groups of an Object
#'
#' @param x The object that contains or is associdated with groups.
#' @param ... Further arguments passed to or from other methods.

#' @return A list or vector of group names. Will return `NULL` if no groups are associated with
#'   the object `x`.
#' @export
#'
#' @family master functions
#'
#'
groups <- function(x, ...) UseMethod("groups")

#' Default Function for Generic `groups`.
#'
#' @param x The object that contains or is associdated with groups.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Will raise an error (`stop`).
#' @export
#'
#' @family master functions
#'
groups.default <- function(x, ...) stop("groups is not implemented for this type")

#' Fetches the Groups defined in the Master Object
#'
#' @param x The master object
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A vector of group names or `NULL` if no groups were defined.
#' @export
#'
#' @family master functions
#'
#' @examples
#'
#' mol <- molecule("Mol 1", MW = 100)
#' m1 <- master("test", mol)
#' m2 <- master("test 2", mol, "my ref", groups = 3)
#' groups(m1) # returns NULL
#' groups(m2) # returns c("Group_1", "Group_2", "Group_3")
#'
groups.master <- function(x, ...) {

  if (!is_master(x))
    stop("Input must be of class master", call. = F)

  return(x$group.names)
}


#' Checks if a Data Structure is a Master Object
#'
#' @param ds The data structure that should be tested.
#'
#' @return True if the data structure is a Master object
#' @export
#'
#' @family master functions
#'
is_master <- function(ds) {
  inherits(ds, "master")
}

#' Adds a New Entry to an Master Object
#'
#' @param pop.file A file or entry descriptor for a population simulation
#'   (can be `NA` or a single string).
#' @param pop.id Id for a population simulation (can be `NA` or a single string).
#'   Default: `pop.file`.
#' @param sim.file A file or entry descriptor for an individual simulation
#'   (can be `NA` or a single string).
#' @param sim.id Id for an individual simulation (can be `NA` or a single string). Default: `pop.file`.
#' @param dose A dose associated with the simulations. Can be `NA` or a single
#'   numeric entry with attached dosing unit. The value will be converted to the dosing unit
#'   associated with the `master.obj`.
#' @param master.obj The master object (must be one valid master object).
#' @param pop.molecules One or more (molecule list) molecules that are associated to
#'   the population simulation. Can be `NA` if `pop.id`/`pop.file` is `NA`.
#' @param sim.molecules One or more (molecule list) molecules that are associated to the
#'   individual simulation. Can be `NA` if `sim.id`/`sim.file` is `NA`.
#' @param obs.ids The associated observed data ids. Can be `NA` or one/more
#'   string-convertibles (e.g. strings or numbers). The Ids are converted to strings.
#' @param groups Group entries for the simulations. Must be `NA` or
#'   a list-like object with string-convertibles or `NA` (if not set for a group).
#' The number of elements cannot be larger than the number of groups of the master object.
#' @param plot.header The name of the plottings header for the simulations.
#'   Can be `NA` or a (empty) single string.
#' @param plot.x.range The plotting range for the x-axis. Can be `NA` (no range).
#'   If not `NA`, units must be associated with the numeric vector of size 2.
#' @param plot.y.range The plotting range for the y-axis. Can be `NA` (no range).
#'   If not `NA`, units must be associated with the numeric vector of size 2.
#'
#' @param entry.list A single or a list of entries of class `masterEntry`. If this parameter
#'   is set all other parameters will be ignored.
#'
#' A note on the plot.X/Y.range options: If only one parameter (with unit) is provided the range
#'   is implicitly set to `c(value, NA)` (no upper limit).
#'   If `NA` is set (one value; no limits) no unit must be set. If only one limit (upper or lower)
#'   should be set the complimentary direction must be set to `NA`.
#'
#' A note on molecules: Molecules must be present in the master object (when master was created).
#'   Unknown molecules will raise an error.
#'
#' If `pop.file` is set `pop.id` must be set as well. If `sim.file` is set `sim.id` must be
#'   set as well. A warning will be displayed if entries with the same `pop.id` or `sim.id` are
#'   already present in the `master.obj`.
#'
#' @return The master object with added entry.
#' @export
#'
#' @family master functions
#'
#' @examples
#'
#' mol_1 <- molecule("Foo", MW = 100)
#' mol_2 <- molecule("Bar", MW = 200)
#' # master with 2 associated molecules
#' m <- master("Test", molecule_list(mol_1, mol_2), groups = 2)
#'
#' m <- add_master_entry(m, pop.file = "Pop", pop.molecules = mol_1) # ok
#' m <- add_master_entry(m, pop.file = "Pop 2", pop.molecules = molecule_list(mol_1, mol_2)) # ok
#'
#' # m <- add_master_entry(m, pop.file = "Pop 2", pop.molecules = NA) # error
#' # error
#' # m <- add_master_entry(m, pop.file = "Pop 2", pop.molecules = molecule("Unknown", MW = 10))
#'
#' # ok only Group_1 is set
#' m <- add_master_entry(m, pop.file = "Pop 5", pop.molecules = mol_1, groups = c("Foo", NA))
#'
#' # ok both groups set
#' m <- add_master_entry(m, pop.file = "Pop 7", pop.molecules = mol_1, groups = c("Foo", "Boo"))
#'
#' # error
#' # m <- add_master_entry(m, pop.name = "Pop 10", pop.molecules = mol_1,
#' #                       groups = c("Foo", "Boo", "Not in master"))
#'
#' m <- add_master_entry(m, plot.x.range = NA) # ok -> no limit
#' m <- add_master_entry(m, plot.x.range = units::as_units(c(1,2), "h")) # ok upper and lower limit
#' m <- add_master_entry(m, plot.x.range = units::as_units(c(NA, 2), "h")) # ok only upper limit
#' m <- add_master_entry(m, plot.x.range = units::as_units(c(1, NA), "h")) # ok only lower limit
#'
#' # m <- add_master_entry(m, plot.x.range = c(1,2)) # error -> no units provided
#'
add_master_entry <- function(master.obj,
                             pop.file = NA,
                             pop.id = pop.file,
                             pop.molecules = NA,
                             sim.file = NA,
                             sim.id = sim.file,
                             sim.molecules = NA,
                             dose = NA,
                             obs.ids = NA,
                             groups = NA,
                             plot.header = NA,
                             plot.x.range = NA,
                             plot.y.range = NA,
                             entry.list = NA
                             ) {

  if (is_masterEntry(entry.list) || is_masterEntry_list(entry.list)) {
    if (is_masterEntry(entry.list)) {
      entry.list <- list(entry.list)
    }

    if (!is_masterEntry_list(entry.list)) {
      stop("entries must be of class masterEntry or a list of masterEntries", call. = F)
    }

    for (x in entry.list) {
      master.obj <- add_master_entry(master.obj,
                                     pop.file = x$Pop_File,
                                     pop.id = x$Pop_Id,
                                     pop.molecules = x$Pop_Mols,
                                     sim.file = x$Sim_File,
                                     sim.id = x$Sim_Id,
                                     sim.molecules = x$Sim_Mols,
                                     dose = x$Dose,
                                     obs.ids = x$Obs_Ids,
                                     groups = x$Groups,
                                     plot.header = x$Plot$Plot_Header,
                                     plot.x.range = x$Plot$X_range,
                                     plot.y.range = x$Plot$Y_range)
    }
    return(master.obj)
  }

  pop.id <- trimws(pop.id)
  sim.id <- trimws(sim.id)
  pop.file <- trimws(pop.file)
  sim.file <- trimws(sim.file)
  obs.ids <- trimws(unlist(obs.ids))
  groups <- trimws(unlist(groups))
  plot.header <- trimws(plot.header)

  if (!is_master(master.obj))
    stop("master.obj must be of class master", call. = F)

  ################################## NAMES ########################################
  if (!is_single_string(pop.id) && !any_na(pop.id))
    stop("pop.id must be a single non-empty string or NA", call. = F)

  if (!is_single_string(pop.file) && !any_na(pop.file))
    stop("pop.file must be a single non-empty string or NA", call. = F)

  if (!is_single_string(sim.id) && !any_na(sim.id))
    stop("sim.id must be a single non-empty string or NA", call. = F)

  if (!is_single_string(sim.file) && !any_na(sim.file))
    stop("sim.file must be a single non-empty string or NA", call. = F)

  # test for both xor set
  if (xor(any_na(pop.id), any_na(pop.file))) {
    stop("pop.id and pop.file must be both NA or both set", call. = F)
  }

  if (xor(any_na(sim.id), any_na(sim.file)))
    stop("sim.id and sim.file must be both NA or both set", call. = F)

  has_pop_entry <- !any_na(pop.id)
  has_sim_entry <- !any_na(sim.id)

  # test for double entries for sim and pop (ids)
  if (!any_na(sim.id) && sim.id %in% master.obj$data$Sim_Id)
    warning(paste("`sim.id` <", sim.id,
                  "> already present in the master object. This is most likely not intended."),
            call. = F)

  if (!any_na(pop.id) && pop.id %in% master.obj$data$Pop_Id)
    warning(paste("`pop.id` <", pop.id,
                  "> already present in the master object. This is most likely not intended."),
            call. = F)

  ################################## MOLS AND DOSE ########################################
  if (is_molecule(pop.molecules))
    pop.molecules <- molecule_list(pop.molecules)
  has_pop_mols <- is_molecule_list(pop.molecules)
  if (xor(has_pop_entry, has_pop_mols))
    stop("pop.molecules and population entries do not match", call. = F)

  if (has_pop_mols && !is_unique_molecule_list(pop.molecules))
    stop("pop.molecules must be unique - double entries found", call. = F)

  if (is_molecule(sim.molecules))
    sim.molecules <- molecule_list(sim.molecules)
  has_sim_mols <- is_molecule_list(sim.molecules)
  if (xor(has_sim_entry, has_sim_mols))
    stop("sim.molecules and individual simulation entries do not match", call. = F)

  if (has_sim_mols && !is_unique_molecule_list(sim.molecules))
    stop("sim.molecules must be unique - double entries found", call. = F)

  if (has_pop_mols && !has_molecules(master.obj$molecules, molecule_ids(pop.molecules)))
    stop("pop.molecules contains unknown molecules", call. = F)

  if (has_sim_mols && !has_molecules(master.obj$molecules, molecule_ids(sim.molecules)))
    stop("sim.molecules contains unknown molecules", call. = F)

  if (!is_single_numeric(dose) && !any_na(dose))
    stop("dose must be NA or single numeric", call. = F)

  if (!any_na(dose) && !has_dose(dose))
    stop("dose must have a valid dosing unit attached", call. = F)

  ################################## GROUPS ########################################
  if (!all(sapply(groups, function(x) is_single_string(x, can.be.empty = F) || any_na(x))))
    stop("groups must be NA or vector of non-empty strings/NA", call. = F)

  if (length(groups) == 1 && any_na(groups))
    groups <- NULL

  if (length(groups) > length(master.obj$group.names))
    stop("number of groups is larger than number of master.obj groups", call. = F)

  entry_groups <- c(groups, rep(NA, length(master.obj$group.names) - length(groups)))

  ################################## HEADER ########################################
  if (!is_single_string(plot.header) && !any_na(plot.header))
    stop("plot.header must be a single non-empty string or NA", call. = F)

  ################################## UNITS ########################################
  if (length(plot.x.range) > 2)
    stop("plot.x.range must be NA or vector of length 2", call. = F)

  has_x_unit <- has_units(plot.x.range)
  x_unit <- if (has_x_unit) units::deparse_unit(plot.x.range) else NA_character_
  plot.x.range <- units::set_units(plot.x.range, NULL)

  # we have to remove the unit first to let sapply work
  if (!all(sapply(plot.x.range, function(x) is_single_numeric(x) || any_na(x))))
    stop("plot.x.range vector elements must be NA or numeric", call. = F)

  has_x_min <- !is.na(plot.x.range[1])
  has_x_max <- !is.na(plot.x.range[2])
  if ((has_x_min || has_x_max) && !has_x_unit)
    stop("plot.x.range must have a unit attached", call. = F)

  x_min <- if (!has_x_min) NA_real_ else plot.x.range[1]
  x_max <- if (!has_x_max) NA_real_ else plot.x.range[2]

  if (length(plot.y.range) > 2)
    stop("plot.y.range must be NA or vector of length 2", call. = F)

  has_y_unit <- has_units(plot.y.range)
  y_unit <- if (has_y_unit) units::deparse_unit(plot.y.range) else NA_character_
  plot.y.range <- units::set_units(plot.y.range, NULL)

  # we have to remove the unit first to let sapply work
  if (!all(sapply(plot.y.range, function(x) is_single_numeric(x) || is.na(x))))
    stop("plot.y.range vector elements must be NA or numeric", call. = F)

  has_y_min <- !is.na(plot.y.range[1])
  has_y_max <- !is.na(plot.y.range[2])
  if ((has_y_min || has_y_max) && !has_y_unit)
    stop("plot.y.range must have a unit attached", call. = F)

  y_min <- if (!has_y_min) NA_real_ else plot.y.range[1]
  y_max <- if (!has_y_max) NA_real_ else plot.y.range[2]

  ################################## OBS ########################################
  obs_has_na <- any_na(obs.ids)
  obs_has_null <- any_null(obs.ids) || length(obs.ids) == 0
  is_string_like <- is_string_list(obs.ids,  string.like = T)
  if (obs_has_null || (obs_has_na && length(obs.ids) > 1) || (!obs_has_na && !is_string_like))
    stop("obs.ids must be NA or singel value or vector of string convertables", call. = F)

  entry <- list()
  entry <- c(entry,
    pop.id,
    pop.file,
    if (has_pop_mols) paste(molecule_ids(pop.molecules), collapse = ",") else NA_character_,
    sim.id,
    sim.file,
    if (has_sim_mols) paste(molecule_ids(sim.molecules), collapse = ",") else NA_character_,
    if (is.na(dose)) NA_real_ else units::set_units(dose, master.obj$dose.unit, mode = "standard"),
    if (obs_has_na) NA_character_ else paste(obs.ids, collapse = ","),
    unlist(entry_groups),
    plot.header,
    x_min,
    x_max,
    x_unit,
    y_min,
    y_max,
    y_unit
  )

  # rbind to data.frame
  names(entry) <- colnames(master.obj$data)
  master.obj$data <- dplyr::bind_rows(master.obj$data, entry)
  return(master.obj)
}

#' Generic Function to  Extract One or More entries (Row Entries) from an Object
#'
#' @param x The object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list of entries.
#'
#' @export
#'
#' @family master functions
#'
entry <- function(x, ...) UseMethod("entry")

#' Default implementation Generic Function to  Extract One or More entries (Row Entries) from an Object
#'
#' @param x The object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Raises an error (`stop`).
#'
#' @export
#'
#' @family master functions
#'
entry.default <- function(x, ...) stop("entry is not implemented for this type")

#' Extracts One or More entries (Row Entries) from a Master Object
#'
#' @param x The master object.
#' @param idx One or more (`vector`) numerical idexes (row numbers) of the master object data.
#' @param ... Further arguments passed to or from other methods.
#'
#' `idx` behaves like `vector[idx]`. E.g. `idx = -1` is allowed and does what the user would
#'   expect.
#'
#' @return A list of entries and deserialized molecules, ids and ranges.
#'  Entries have class `masterEntry`. If only one object is fetched only one `masterEntry` object
#'  is returned.
#'
#' @export
#'
#' @family master functions
#'
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' mol <- molecule("Mol 1", MW = 100)
#' m1 <- master("test", mol)
#' m1 <- add_master_entry(m1, "Pop 1", pop.molecules = mol)
#' m1 <- add_master_entry(m1, "Pop 2", pop.molecules = mol)
#'
#' m_filtered <- entry(m1, 1) # one entry (not a list of entries)
#'
entry.master <- function(x, idx, ...) {

  if (!is_master(x))
    stop("master.obj must be of class master", call. = F)

  if (!idx %>% purrr::every(is_single_numeric)) {
    stop("idx must be numeric", call. = F)
  }

  tmp_fn <- function(row) {
    res <- list( Pop_Id = row$Pop_Id,
            Pop_File = row$Pop_File,
            Pop_Mols = if (any_na(row$Pop_Mols))
              NA
            else
              molecule_from_ids(x$molecules, query.ids = unlist(strsplit(row$Pop_Mols, ","))),
            Sim_Id = row$Sim_Id,
            Sim_File = row$Sim_File,
            Sim_Mols = if (any_na(row$Sim_Mols))
              NA
            else
              molecule_from_ids(x$molecules, query.ids = unlist(strsplit(row$Sim_Mols, ","))),
            Dose = if (any_na(row$Dose)) NA else units::as_units(row$Dose, x$dose.unit),
            Obs_Ids = if (any_na(row$Obs_Ids)) NA else unlist(strsplit(row$Obs_Ids, ",")),
            Plot = list(
              Plot_Header = row$Plot_Header,
              X_range = to_range(row$Plot_X_Min, row$Plot_X_Max, row$Plot_X_Unit),
              Y_range = to_range(row$Plot_Y_Min, row$Plot_Y_Max, row$Plot_Y_Unit)
            ),
            Groups = if (is.null(x$group.names)) NULL else row[[x$group.names]])

    class(res) <- append("masterEntry", class(res))
    return(res)
  }

  dat <- x$data[idx,]

  result <- list()
  for (i in 1:nrow(dat)) {
    result <- result %>% append(list(tmp_fn(dat[i,])))
  }

  if (length(idx) == 1)
    result <- result[[1]]

  return(result)
}

#' Checks if a Data Structure is a MasterEntry Object
#'
#' @param x The data structure that should be tested.
#'
#' @return True if the data structure is a MasterEntry Object
#' @export
#'
#' @family master functions
#'
is_masterEntry <- function(x) {
  inherits(x, "masterEntry")
}

#' Checks if a Data Structure is a list of MasterEntry Objects
#'
#' @param list The data structure that should be tested.
#'
#' @return `TRUE` if the data structure is a MasterEntry list. For single MasterEntry
#'   this function returns `FALSE`.
#' @export
#'
#' @family master functions
#'
is_masterEntry_list <- function(list) {

  if (length(list) < 1)
    return(F)

  all(sapply(list, is_masterEntry))
}


#' Extracts One or More entries (Row Entries) from a Master Object
#'
#' @param x The master object.
#' @param ... One or more (`vector`) numerical idexes (row numbers) of the master object data.
#'
#' `...` behaves like `vector[idx]`. E.g. `idx = -1` is allowed and does what the user would
#'   expect.
#'
#' @return A list of entries and deserialized molecules, ids and ranges.
#'  Entries have class `masterEntry`. If only one object is fetched only one `masterEntry` object
#'  is returned.
#'
#' @export
#'
#' @family master functions
#'
#' @examples
#'
#' mol <- molecule("Mol 1", MW = 100)
#' m1 <- master("test", mol)
#' m1 <- add_master_entry(m1, "Pop 1", pop.molecules = mol)
#' m1 <- add_master_entry(m1, "Pop 2", pop.molecules = mol)
#'
#' m_filtered <- m1[1] # one entry (not a list of entries)
#' m_filtered_2 <- m1[1,2] # two entries
#'
`[.master` <- function(x, ...) {
  l <- list(...)
  entry(x, unlist(l))
}


# TODO: Tests and examples
#' Binds (Row Binds) two Master Objects
#'
#' @param master The master object that meta data is preserved.
#' @param slave The slave master object that is bound to the master object.
#'   The number of groups of the slave object must be at most the number of groups of the master
#'   object. Missing group entries will be filled with `NA`.
#'
#' @return A bound (combined) master object.
#' @export
#'
#' @family master functions
#'
#' @examples
#'
#' mol_1 <- molecule("Foo", MW = 100)
#' mol_2 <- molecule("Bar", MW = 200)
#' a <- master("master", mol_1, groups = 2)
#' b <- master("slave", mol_2, groups = 1)
#' bind_master(a, b)
#' # bind_master(b, a) # fails
#'
bind_master <- function(master, slave) {

  if (!is_master(master))
    stop("master must be of class master")

  if (!is_master(slave))
    stop("slave must be of class master")

  if (length(master$group.names) < length(slave$group.names))
    stop("number of master groups must be at least the number of slave groups")


  if (ncol(master$data) != ncol(slave$data)) {
    n.miss <- ncol(master$data) - ncol(slave$data)
    fill <- as.data.frame(matrix(NA, nrow(slave$data), n.miss))

    s.stop.1 <- 8 + length(slave$group.names)
    s.stop.2 <- ncol(slave$data)
    slave$data <- cbind(slave$data[1:s.stop.1], fill, slave$data[(s.stop.1 + ncol(fill)):s.stop.2])
  }

  # just be be sure and for the unequal grouping situation
  colnames(slave$data) <- colnames(master$data)

  if (master$dose.unit != slave$dose.unit) {
    slave$data$Dose <- convert_values(slave$data$Dose,
                                      from = slave$dose.unit,
                                      to = master$dose.unit,
                                      attach.unit = F)
  }
  master$data <- rbind(master$data, slave$data)

  # test for double entries for sim and pop (ids)
  if (any(duplicated(master$data$Sim_Id, incomparables = NA)))
    warning("Found duplicated Sim_Ids in the master object. This is most likely not intended.",
            call. = F)

  if (any(duplicated(master$data$Pop_Id, incomparables = NA)))
    warning("Found duplicated Pop_Ids in the master object. This is most likely not intended.",
            call. = F)

  # molecule merge
  all_mols <- c(master$molecules, slave$molecules)
  all_mol_ids <- molecule_ids(all_mols)
  master$molecules <- molecule_from_ids(all_mols, unique(all_mol_ids))
  if (is_molecule(master$molecules))
    master$molecules <- molecule_list(master$molecule)

  return(master)
}

# TODO: Tests and examples
#' Filters a Master Object
#'
#' @param master.obj The master object that should be filtered.
#' @param entries The entries (rows) that should be returned.
#'   `NA` if this option is not used. Expects a numeric vector.
#' @param valid.pop If valid (non-`NA`) population profile entries should
#'   be returned. `NA` if this option is not used. Expects a logical.
#' @param valid.sim If valid (non-`NA`) individual simulation profile entries shoudl be returned.
#'   `NA` if this option is not used. Expects a logical.
#' @param group.fiter If not NA is provided a group filter is applied
#'   (see details for further information).
#' @param inv.group.filter Logical value to indicate if the group filter should be inverted.
#'
#' Notes on the group filter: If `NA` is set, no group filtering is performed.
#'   To use a group filter the function expects a named list with string vector
#'   values for the matching. The name of a list entry must be one of the group names
#'   of the master object. All group filter are combined with a logical AND.
#'   The values of a list entry will match exactly for the corresponding group entries
#'   of the master object and are combined with a logical OR.
#'
#' All filters are combined with a logical AND.
#'
#' @return The filtered master object.
#' @export
#'
#' @family master functions
#'
#' @examples
#'
#' mol_1 <- molecule("Foo", MW = 100)
#' mol_2 <- molecule("Bar", MW = 200)
#' mol_3 <- molecule("Baz", MW = 300)
#' m <- master("master", molecules = molecule_list(mol_1, mol_2, mol_3), groups = 3)
#' m <- add_master_entry(m, "Pop", pop.molecules = mol_1, groups = c(NA, "AB", "B"))
#' m <- add_master_entry(m, sim.file = "Sim", sim.molecules = mol_2, groups = c("C", "B", "B"))
#' m <- add_master_entry(m, sim.file = "Sim 2", sim.molecules = mol_2, groups = c(NA, "AB", NA))
#' m <- add_master_entry(m, sim.file = "Sim 3",
#'                       sim.molecules = molecule_list(mol_1, mol_2), groups = c(NA, "AB", NA))
#'
#' filter_master(m, entries = c(1,4)) # filter for rows
#' # filter inverse for groups
#' filter_master(m, group.fiter = list(Group_2 = "B", Group_3 = "B"), inv.group.filter = TRUE)
#' filter_master(m, valid.pop = TRUE) # only valid population profiles
#' filter_master(m, valid.sim = TRUE) # only valid individual simulation profiles
#'
filter_master <- function(master.obj,
                          entries = NA,
                          valid.pop = NA,
                          valid.sim = NA,
                          group.fiter = NA,
                          inv.group.filter = F) {

  if (!is_master(master.obj))
    stop("master.obj must be of class master")

  if (!is.numeric(entries) && !single_na(entries))
    stop("entries must be a single NA or numeric")

  if (!is_single_logical(valid.pop) && !single_na(valid.pop))
    stop("valid.pop must be NA or logical")

  if (!is_single_logical(valid.sim) && !single_na(valid.sim))
    stop("valid.sim must be NA or logical")

  if (!is.list(group.fiter) && !single_na(group.fiter))
    stop("group.filter must be NA or named list with no NA values")

  if (is.list(group.fiter)) {
    g_names <- names(group.fiter)
    if (!all(g_names %in% master.obj$group.names))
      stop("Unknown group names defined")
  }

  if (!is_single_logical(inv.group.filter))
    stop("inv.group.filter must be logical")

  data <- master.obj$data

  # rows
  if (is.numeric(entries)) {
    data <- data[entries,]
  }

  # valid pops
  if (!is.na(valid.pop)) {
    take_rows <- if (valid.pop) !is.na(data$Pop_Id) else is.na(data$Pop_Id)
    data <- data[take_rows,]
  }

  # valid sims
  if (!is.na(valid.sim)) {
    take_rows <- if (valid.sim) !is.na(data$Sim_Id) else is.na(data$Sim_Id)
    data <- data[take_rows,]
  }

  # group filter
  if (is.list(group.fiter)) {
    g_names <- names(group.fiter)
    for (group in g_names) {
      filters <- group.fiter[[group]]
      take_rows <- if (inv.group.filter) !(data[[group]] %in% filters) else data[[group]] %in% filters
      data <- data[take_rows,]
    }
  }

  master.obj$data <- data
  return(master.obj)
}


#' Prints Information about a Master Object
#'
#' @param x The Master object.
#' @param details True if all details should be printed, else False.
#' @param max Maximum number of data rows to print. If `max = NA` (default)
#'   the number of entries that are printed
#'  is `getOption("max.print")/ncol(data)`
#'  for detailed print and the number of rows printed is 10 for non-detailed print.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisible copy of the object `x`.
#'
#' @family master functions
#'
#' @export
#'
#' @examples
#' a <- master("Master 1", molecules = molecule("Test", MW = 100))
#' a # prints a without details
#' print(a, details = TRUE) # prints a with details
#' print(a, max = 5) # prints a without details and maximum of 5 data rows
#'
print.master <- function(x, details = F, max = NA, ...) {

  if (!is_master(x))
    stop("Input must be of class master")

  if (!is_single_logical(details))
    stop("details must be logical")

  if (!is_single_numeric(max) && !is.na(max) && max < 0)
    stop("max must be a positive number or NA")

  if (is.na(max))
    max <- if (details) getOption("max.print")/ncol(x$data) else 10

  cat(paste(x$id, "(Master)\n"))
  cat(paste("Reference:\t", x$reference, "\n"))
  cat(paste("Molecules:\t", paste(molecule_ids(x$molecules), collapse = ", "), "\n"))
  cat(paste("Dose unit:\t", x$dose.unit , "\n"))
  cat(paste("Groups:\t\t", paste(x$group.names, collapse = ", "), "\n"))
  cat(paste("Entries:\t", nrow(x$data), "\n"))

  if (details) {
    print(x$data, row.names = F, n = max)
  }
  else {
    dat <- x$data[,1:(8 + length(x$group.names))]
    print(dat, row.names = F, n = max)
  }

  invisible(x)
}
