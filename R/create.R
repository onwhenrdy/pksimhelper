
# creates a molecule
molecule <- function(name,
                     display.name = NA,
                     id = NA,
                     file.name.match = NA,
                     add.file.matcher = c(),
                     pubchem.id = NA,
                     MW = NA,
                     is.fraction = F,
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
