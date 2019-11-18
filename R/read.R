
# extract time column information and unit
.id.time.col <- function(df, match = "time [") {
  col.names <- tolower(colnames(df))
  col.nr <- which(grepl(match, col.names, fixed = TRUE))

  if (length(col.nr) != 1) {
    stop("Could not identify time column")
  }

  result <- list(id = "Time", col = col.nr, unit = .extract.unit(col.names[col.nr]))

  return(result)
}


# extract molecule column information and unit
.id.molecule.cols <- function(df, match.tag, add.match.tags = c(), is.fraction = F, fixed.unit = NA,
                              silent = F) {
  col.names <- tolower(colnames(df))
  match.tag <- tolower(match.tag)
  if (length(add.match.tags) > 0) {
    add.match.tags <- tolower(add.match.tags)
  }

  grep.tag <- match.tag

  if (length(add.match.tags > 0)) {
    grep.tag <- c(grep.tag, add.match.tags)
    grep.tag <- sapply(grep.tag, function(x) paste("(?=.*", x, ")", sep = ""), simplify = T, USE.NAMES = F)
    grep.tag <- paste(grep.tag, collapse = "")
  }

  col.nrs <- c(grep(grep.tag, col.names, perl = TRUE))
  n.cols <- length(col.nrs)

  if (n.cols < 1) {
    if (silent)
      return(NULL)

    stop(paste("Could not identify value column for match tag <", match.tag, ">"))
  }

  unit <- NA
  if (!is.na(fixed.unit)) {
    unit <- fixed.unit
  } else if (is.fraction) {
    unit <- units::as_units(units::unitless)
  } else {
    unit <- .extract.unit(col.names[col.nrs[1]])
  }

  result <- list(
    id = match.tag,
    cols = col.nrs,
    unit = unit
  )

  return(result)
}


.id.col <- function(df, match, fixed = TRUE) {
  col.names <- trimws(tolower(colnames(df)))
  match <- tolower(match)

  col.nr <- which(grepl(match, col.names, fixed = fixed))

  if (identical(col.nr, integer(0))) {
    return(NA)
  }

  return(col.nr)
}


.match.mols.from.string <- function(ids, molecules) {
  mols <- list()
  ids <- trimws(tolower(ids))
  for (id in ids) {
    match <- FALSE
    for (mol in molecules) {
      if (id == tolower(mol$id)) {
        match <- TRUE
        mols <- append(mols, list(mol))
        break
      }
    }

    if (!match) {
      stop(paste("Could not find molecule with id <", id, ">"))
    }
  }

  return(mols)
}

groups.master <- function(master.data) {
  group_idx <- which(grepl("group", tolower(colnames(master.data$data))))
  has_groups <- length(group_idx) > 0
  if (has_groups) {
    return(colnames(master.data$data)[group_idx])
  }

  return(NULL)
}


read.pop.profiles.from.master <- function(master.data,
                                          molecules = list(),
                                          files.dir = ".",
                                          pop.file.ext = "",
                                          pop.file.format = c("auto", "xsl", "csv"),
                                          pop.file.csv.sep = ",",
                                          pop.file.csv.dec = ".",
                                          pop.file.xls.sheet = 1) {
  if (!is.master(master.data)) {
    stop("Input must be of class master")
  }

  data <- master.data$data
  lines <- which(!is.na(data$pop))

  group_names <- groups.master(master.data)

  results <- list()
  for (i in lines) {
    pop <- trimws(data$pop[i])
    mols <- trimws(data$pop.mol[i])
    obs.ids <- if (is.na(data$obs[i])) NA else trimws(unlist(strsplit(data$obs[i], ",")))
    mols <- .match.mols.from.string(unlist(strsplit(mols, ",", fixed = TRUE)),
      molecules = molecules
    )

    fixed.pop <- gsub("/", "_", pop)
    file <- paste0(fixed.pop, pop.file.ext)
    file <- file.path(files.dir, file)

    profiles <- read.pop.profiles(file,
      molecules = mols,
      format = pop.file.format,
      csv.sep = pop.file.csv.sep,
      csv.dec = pop.file.csv.dec,
      xls.sheet = pop.file.xls.sheet
    )

    plot.infos <- list(
      x.min = data$x.min[i],
      x.max = data$x.max[i],
      y.min = data$y.min[i],
      y.max = data$y.max[i],
      y.min = data$y.min[i],
      y.max = data$y.max[i],
      y.min.log = data$y.min.log[i],
      y.max.log = data$y.max.log[i],
      x.offset = data$x.offset[i],
      main = data$main[i]
    )

    groups <- NULL
    if (length(group_names) > 0)
      groups <- as.character(data[i, group_names])

    pop.sim <- list(
      id = pop, obs.ids = obs.ids, profiles = profiles,
      groups = groups,
      plot.infos = plot.infos
    )
    class(pop.sim) <- "MatchedProfiles"
    results <- append(results, list(pop.sim))
  }

  return(results)
}

is.matched.profiles <- function(x) inherits(x, "MatchedProfiles")


# Read study list from file
read.pop.profiles <- function(file,
                              molecules = list(),
                              format = c("auto", "xsl", "csv"),
                              add.file.ext = "",
                              csv.sep = ",",
                              csv.dec = ".",
                              xls.sheet = 1) {
  file <- paste0(file, add.file.ext)
  base.name <- basename(file)
  message(paste("Reading profile from file <", base.name, ">"))

  if (!.does.file.exist(file)) {
    stop(paste("File <", file, "> does not exist"))
  }

  if (!.is.file.readable(file)) {
    stop(paste("File <", file, "> is not readable"))
  }

  format <- match.arg(format)

  if (format == "auto") {
    format <- .identify.file.ext(file)
    if (is.na(format)) {
      stop(paste("File <", file, "> has unknown file extension. Please specify the format."))
    }
  }

  if (format == "tsv") {
    format <- "csv"
  }

  df <- NA
  if (format == "xls") {
    excel.sheet <- suppressMessages(readxl::read_excel(file,
      sheet = xls.sheet,
      col_names = TRUE
    ))
    df <- as.data.frame(excel.sheet)
  }

  if (format == "csv") {
    df <- .read.csv(
      file = file, header = TRUE, sep = csv.sep, dec = csv.dec,
      encoding = "UTF-8"
    )
  }

  # extract time, individual and molecule columns
  time.info <- .id.time.col(df)
  mol.infos <- list()
  for (mol in molecules) {
    match <- .id.molecule.cols(df, mol$file.name.match,
      mol$add.file.matcher,
      is.fraction = mol$is.fraction, fixed.unit = mol$fixed.unit
    )

    if (length(match$cols) > 1) {
      stop(paste("Ambiguous column found for tag <", mol$file.name.match, ">"))
    }

    mol.infos <- c(mol.infos, list(match))
  }

  id.col <- .id.col(df, "IndividualId")
  if (is.na(id.col) || length(id.col) > 1) {
    stop("Could not identify individual id column")
  }

  # extract time data and id data
  times <- sort(unique(df[, time.info$col]))
  # prepare results
  results <- list()
  for (i in 1:length(molecules)) {
    mol <- molecules[[i]]
    entry <- list(
      molecule = mol,
      id = "",
      reference = "",
      group = NA,
      time.unit = time.info$unit,
      value.unit = mol.infos[[i]]$unit,
      data = data.frame(Time = times),
      type = "population",
      data.type = "individual",
      origin = "sim"
    )
    class(entry) <- "profile"
    results <- c(results, list(entry))
  }
  # for every id
  ids <- sort(unique(df[, id.col]))
  for (id in ids) {
    id.df <- df[df[id.col] == id, ]
    id.df <- id.df[order(id.df[, time.info$col]), ]
    for (i in 1:length(molecules)) {
      results[[i]]$data <- cbind(results[[i]]$data, id.df[, mol.infos[[i]]$cols])
    }

    # check for duplicated times
    tmp.data <- results[[i]]$data
    dup.times <- duplicated(tmp.data$Time)
    if (any(dup.times)) {
      message(paste("Deleted duplicated time entries for file <", base.name, ">"))
      results[[i]]$data <- tmp.data[!dup.times, ]
    }
  }

  for (i in 1:length(molecules)) {
    colnames(results[[i]]$data) <- c("Time", paste("ID", ids))
  }

  return(results)
}

is.profile <- function(x) inherits(x, "profile")


read.obs.profiles <- function(file,
                              obs.sheet = 1,
                              reference.sheet = 2,
                              molecules = list(),
                              id.filter = NA) {
  if (length(id.filter) > 1 || !is.na(id.filter)) {
    id.filter <- unlist(lapply(id.filter, toString))
  }

  base.name <- basename(file)
  if (!.does.file.exist(file)) {
    stop(paste("File <", base.name, "> does not exist"))
  }

  if (!.is.file.readable(file)) {
    stop(paste("File <", base.name, "> is not readable"))
  }

  if (!.does.sheet.exist(file, obs.sheet)) {
    stop(paste("Observed data sheet in file <", base.name, "> does not exist"))
  }

  if (!.does.sheet.exist(file, reference.sheet)) {
    stop(paste("Observed data sheet in file <", base.name, "> does not exist"))
  }


  # reference sheet
  # we text parsing because auto-guess will freak out on mixed column types
  ref.sheet <- suppressMessages(readxl::read_excel(file,
    sheet = reference.sheet,
    col_names = TRUE,
    col_types = "text"
  ))

  ref.sheet <- as.data.frame(ref.sheet)
  id.col <- .id.col(ref.sheet, "ID")
  if (is.na(id.col)) {
    stop(paste("File <", base.name, ">: Could not identify ID column in reference sheet"))
  }

  mol.col <- .id.col(ref.sheet, "Compound")
  if (is.na(mol.col)) {
    stop(paste("File <", base.name, ">: Could not identify Compound column in reference sheet"))
  }

  ref.col <- .id.col(ref.sheet, "Reference")
  if (is.na(ref.col)) {
    stop(paste("File <", base.name, ">: Could not identify Reference column in reference sheet"))
  }

  group.col <- .id.col(ref.sheet, "Group")
  if (is.na(group.col[1])) {
    stop(paste("File <", base.name, ">: Could not identify Group column in reference sheet"))
  }

    cite.col <- .id.col(ref.sheet, "Citekey")
  if (is.na(cite.col)) {
    stop(paste("File <", base.name, ">: Could not identify Citekey column in reference sheet"))
  }

  dose.col <- .id.col(ref.sheet, "^Dose", fixed = FALSE)
  if (is.na(dose.col)) {
    stop(paste("File <", base.name, ">: Could not identify Dose column in reference sheet"))
  }

  dunit.col <- .id.col(ref.sheet, "Unit Dose")
  if (is.na(dunit.col)) {
    stop(paste("File <", base.name, ">: Could not identify Unit Dose column in reference sheet"))
  }

  route.col <- .id.col(ref.sheet, "Route")
  if (is.na(route.col)) {
    stop(paste("File <", base.name, ">: Could not identify Unit Route column in reference sheet"))
  }
  ref.sheet <- ref.sheet[c(id.col, mol.col, ref.col, group.col, cite.col, dose.col, dunit.col, route.col)]
  colnames(ref.sheet) <- c("ID", "MOL", "REF", "GROUP", "GROUP2", "GROUP3", "CKEY", "DOSE", "DOSEUNIT", "ROUTE")
  ref.sheet$ID <- trimws(ref.sheet$ID)
  ref.sheet$GROUP <- trimws(ref.sheet$GROUP)
  ref.sheet$GROUP2 <- trimws(ref.sheet$GROUP2)
  ref.sheet$GROUP3 <- trimws(ref.sheet$GROUP3)
  id.mol.list <- list()
  for (i in 1:length(ref.sheet$MOL)) {
    mol <- ref.sheet$MOL[i]
    tmp.mol <- .find.molecule.from.id(molecules, mol)
    if (length(tmp.mol) < 2 && is.na(tmp.mol)) {
      stop(paste("File <", base.name, ">: Could not find molecule with id <", mol, "> from the reference sheet"))
    }
    id.mol.list <- c(id.mol.list, list(tmp.mol))
  }

  # checks
  if (length(ref.sheet$ID) != length(unique(ref.sheet$ID))) {
    stop(paste("File <", base.name, ">: Found non-unique IDs in the reference sheet:", ref.sheet$ID[which(duplicated(ref.sheet$ID))]))
  }

  # observed data sheet
  # we text parsing because auto-guess will freak out on mixed column types
  obs.sheet <- suppressMessages(readxl::read_excel(file,
    sheet = obs.sheet,
    col_names = TRUE,
    col_types = "text"
  ))
  df <- as.data.frame(obs.sheet)
  id.col <- .id.col(df, "ID")
  if (is.na(id.col)) {
    stop(paste("File <", base.name, ">: Could not identify ID column for observed sheet"))
  }

  time.col <- .id.time.col(df)
  value.col <- .id.molecule.cols(df, "Value")
  error.col <- .id.molecule.cols(df, "SD")

  df[, time.col$col] <- sapply(df[, time.col$col], as.numeric)
  df[, value.col$cols] <- sapply(df[, value.col$cols], as.numeric)
  df[, error.col$cols] <- sapply(df[, error.col$col], as.numeric)

  ids <- unique(df[, id.col])
  if (length(id.filter) > 1 || !is.na(id.filter)) {
    ids <- id.filter
  }

  results <- list()
  for (id in ids) {
    data <- df[df[, id.col] == id, ]
    data <- data[c(time.col$col, value.col$cols, error.col$cols)]
    data <- cbind(data, Max = data[, 2] + data[, 3])
    data[, 3] <- data[, 2] - data[, 3]
    colnames(data) <- c("Time", "Avg", "Min", "Max")
    i <- which(ref.sheet$ID == id)
    if (length(i) != 1) {
      warning(paste("File <", base.name, ">: Found observed ID <", id,
                    "> that has no match in the reference data sheet. Data will be skipped"))
      next
    }

    # check for duplicated times
    dup.times <- duplicated(data$Time)
    if (any(dup.times)) {
      message(paste(
        "Deleted duplicated time entries for file <",
        base.name, "> and molecule id <", id.mol.list[[i]]$id, ">"
      ))
      data <- data[!dup.times, ]
    }

    entry <- list(
      molecule = id.mol.list[[i]],
      reference = gsub("\r\n", "\n", ref.sheet$REF[i]),
      citekey = trimws(ref.sheet$CKEY[i]),
      dose = trimws(ref.sheet$DOSE[i]),
      dose.unit = trimws(ref.sheet$DOSEUNIT[i]),
      route = trimws(ref.sheet$ROUTE[i]),
      group = ref.sheet$GROUP[i],
      group2 = ref.sheet$GROUP2[i],
      group3 = ref.sheet$GROUP3[i],
      id = id,
      time.unit = time.col$unit,
      value.unit = value.col$unit,
      data = data,
      type = "individual",
      data.type = "mean",
      origin = "obs"
    )

    colnames(ref.sheet) <- c("ID", "MOL", "REF", "GROUP", "GROUP2", "GROUP3", "CKEY", "DOSE", "DOSEUNIT", "ROUTE")


    class(entry) <- "profile"

    results <- c(results, list(entry))
  }

  message(paste("Parsed", length(ids), "oberved profiles"))
  return(results)
}


read.master.file <- function(master.file,
                             molecules = list(),
                             format = c("auto", "xsl", "csv"),
                             csv.sep = ",",
                             csv.dec = ".",
                             xls.sheet = 1) {
  base.name <- basename(master.file)

  if (!.does.file.exist(master.file)) {
    stop(paste("File <", base.name, "> does not exist"))
  }

  if (!.is.file.readable(master.file)) {
    stop(paste("File <", base.name, "> is not readable"))
  }

  format <- match.arg(format)

  if (format == "auto") {
    format <- .identify.file.ext(master.file)
    if (is.na(format)) {
      stop(paste("File <", base.name, "> has unknown file extension. Please specify the format."))
    }
  }

  if (format == "tsv") {
    format <- "csv"
  }

  df <- NA
  if (format == "xls") {
    suppressMessages(excel.sheet <- readxl::read_excel(master.file,
      sheet = xls.sheet,
      col_names = TRUE
    ))
    df <- as.data.frame(excel.sheet)
  }

  if (format == "csv") {
    df <- .read.csv(
      file = master.file, header = TRUE, sep = csv.sep, dec = csv.dec,
      encoding = "UTF-8"
    )
  }

  c.names <- tolower(colnames(df))
  # column identifier
  pop.id <- which(grepl("pop_name", c.names))
  sim.id <- which(grepl("sim_name", c.names))
  obs.id <- which(grepl("obs", c.names))
  pop.molecules.id <- which(grepl("pop_mol", c.names))
  sim.molecules.id <- which(grepl("sim_mol", c.names))

  if (length(pop.id) == 0 || length(pop.id) > 1) {
    stop(paste("File <", base.name, "> has ambiguous or unknown population column."))
  }

  if (length(sim.id) == 0 || length(sim.id) > 1) {
    stop(paste("File <", base.name, "> has ambiguous or unknown simulation column."))
  }

  if (length(obs.id) == 0 || length(obs.id) > 1) {
    stop(paste("File <", base.name, "> has ambiguous or unknown observed column."))
  }

  if (length(pop.molecules.id) == 0 || length(pop.molecules.id) > 1) {
    stop(paste("File <", base.name, "> has ambiguous or unknown population molecule column."))
  }

  if (length(sim.molecules.id) == 0 || length(sim.molecules.id) > 1) {
    stop(paste("File <", base.name, "> has ambiguous or unknown simulation molecule column."))
  }

  # optionals
  xmin.id <- which(grepl("x_min", c.names))
  if (length(xmin.id) == 1) {
    units(df[, xmin.id]) <- .extract.unit(c.names[xmin.id])
    c.names[xmin.id] <- "x.min"
  }
  xmax.id <- which(grepl("x_max", c.names))
  if (length(xmax.id) == 1) {
    units(df[, xmax.id]) <- .extract.unit(c.names[xmax.id])
    c.names[xmax.id] <- "x.max"
  }

  ymin.id <- which(grepl("^y_min", c.names))
  if (length(ymin.id) == 1) {
    units(df[, ymin.id]) <- .extract.unit(c.names[ymin.id])
    c.names[ymin.id] <- "y.min"
  }
  ymax.id <- which(grepl("^y_max", c.names))
  if (length(ymax.id) == 1) {
    units(df[, ymax.id]) <- .extract.unit(c.names[ymax.id])
    c.names[ymax.id] <- "y.max"
  }

  log_y_min.id <- which(grepl("^log_y_min", c.names))
  if (length(log_y_min.id) == 1) {
    units(df[, log_y_min.id]) <- .extract.unit(c.names[log_y_min.id])
    c.names[log_y_min.id] <- "y.min.log"
  }
  log_ymax.id <- which(grepl("^log_y_max", c.names))
  if (length(log_ymax.id) == 1) {
    units(df[, log_ymax.id]) <- .extract.unit(c.names[log_ymax.id])
    c.names[log_ymax.id] <- "y.max.log"
  }

  offset.id <- which(grepl("offset", c.names))
  if (length(offset.id) == 1) {
    units(df[, offset.id]) <- .extract.unit(c.names[offset.id])
    c.names[offset.id] <- "x.offset"
  }

  headline.id <- which(grepl("headline", c.names))
  c.names[headline.id] <- "main"

  c.names[pop.id] <- "pop"
  c.names[sim.id] <- "sim"
  c.names[obs.id] <- "obs"
  c.names[pop.molecules.id] <- "pop.mol"
  c.names[sim.molecules.id] <- "sim.mol"
  colnames(df) <- c.names

  group_cols <- which(grepl("group", c.names))
  group_cols <- c.names[group_cols]

  # drop unknown colums
  known <- c(
    "pop", "sim", "obs", "pop.mol", "sim.mol",
    "x.min", "x.max",
    "y.min", "y.max",
    "y.min.log", "y.max.log",
    "x.offset", "main", group_cols)
  df <- df[, colnames(df) %in% known]

  if ("obs" %in% colnames(df))
    df$obs <- as.character(df$obs)

  # test for missing molecules
  for (i in 1:nrow(df)) {
    mols <- df$pop.mol
    mol.strs <- unlist(strsplit(mols, ",", fixed = TRUE))
    if (length(mol.strs) == 0) {
      stop(paste("File <", base.name, "> in line <", i, ">: no population molecule definition found."))
    }

    .match.mols.from.string(mol.strs, molecules)

    mols <- df$sim.mol
    mol.strs <- unlist(strsplit(mols, ",", fixed = TRUE))
    if (length(mol.strs) == 0) {
      stop(paste("File <", base.name, "> in line <", i, ">: no simulation molecule definition found."))
    }

    .match.mols.from.string(mol.strs, molecules)
  }

  result <- list(id = base.name, data = df)
  class(result) <- "master"
  return(result)
}

is.master <- function(x) inherits(x, "master")


# read all sheets or csv data into a list of dataframes
.read.all.sheet <- function(files,
                            folder = ".",
                            format = c("auto", "xsl", "csv"),
                            csv.sep = ",",
                            csv.dec = ".") {
  format <- match.arg(format)

  results <- list()
  for (file in files) {
    base.name <- basename(file)
    file <- file.path(folder, file)
    message(paste("Reading file <", base.name, ">"))

    if (!.does.file.exist(file)) {
      stop(paste("File <", base.name, "> does not exist"))
    }

    if (!.is.file.readable(file)) {
      stop(paste("File <", base.name, "> is not readable"))
    }

    f.format <- format
    if (f.format == "auto") {
      f.format <- .identify.file.ext(file)
    }

    if (is.na(f.format)) {
      stop(paste("File <", base.name, "> : Could not detect file format"))
    }

    if (f.format == "tsv") {
      f.format <- "csv"
    }

    df <- NA
    if (f.format == "csv") {
      df <- .read.csv(
        file = file, header = TRUE, sep = csv.sep, dec = csv.dec,
        encoding = "UTF-8"
      )

      attr(df, "file.name") <- base.name
      attr(df, "sheet.idx") <- 1
      attr(df, "sheet.name") <- NA
      results <- append(results, df)
    } else {
      sheet.names <- .sheets(file)
      sheets <- length(sheet.names)
      for (i in 1:sheets) {
        suppressMessages(excel.sheet <- readxl::read_excel(file,
          sheet = i,
          col_names = TRUE
        ))

        df <- as.data.frame(excel.sheet)
        attr(df, "file.name") <- base.name
        attr(df, "sheet.idx") <- i
        attr(df, "sheet.name") <- sheet.names[i]

        results <- append(results, list(df))
      }
    }
  }

  return(results)
}


.get.sheet.info.str <- function(df) {
  file.name <- attributes(df)["file.name"]
  sheet.idx <- attributes(df)["sheet.idx"]
  sheet.name <- attributes(df)["sheet.name"]
  return(paste(
    "  File: ", file.name,
    "\n  Sheet Nr:", sheet.idx,
    "\n  Sheet Name:", sheet.name
  ))
}

read.sim.profiles.from.master <- function(master.data,
                                          sim.files,
                                          molecules = list(),
                                          files.dir = ".",
                                          file.format = c("auto", "xsl", "csv"),
                                          file.csv.sep = ",",
                                          file.csv.dec = ".",
                                          action.on.multimatch = c("stop", "warning", "message", "silent"),
                                          multifile.enties = F
                                          ) {
  if (!is.master(master.data)) {
    stop("Input must be of class master")
  }


  action.on.multimatch <- match.arg(action.on.multimatch)
  multi.action <- stop
  if (action.on.multimatch == "warning") {
    multi.action <- warning
  } else if (action.on.multimatch == "message") {
    multi.action <- message
  } else if (action.on.multimatch == "silent") {
    multi.action <- function(...) invisible(NULL)
  } # noop


  results <- list()
  sim.data <- .read.all.sheet(sim.files,
    folder = files.dir, format = file.format,
    csv.sep = file.csv.sep, csv.dec = file.csv.dec
  )

  master <- master.data$data
  for (i in 1:nrow(master)) {
    id <- master$sim[i]
    if (is.na(id)) {
      next
    }

    message(paste("Lookup id <", id, ">"))

    # match molecules
    mols <- trimws(master$sim.mol[i])
    obs.ids <- if (is.na(master$obs[i])) NA else trimws(unlist(strsplit(master$obs[i], ",")))
    mols <- .match.mols.from.string(unlist(strsplit(mols, ",", fixed = TRUE)),
      molecules = molecules
    )

    # find entry in sim.data
    sheet_entry <- list()
    prev_matches <- NA
    for (sheet in sim.data) {
      matches <- .id.col(sheet, paste0(id, "|"))
      if (length(matches) > 1 || !is.na(matches)) {
        # found it and check for double entry
        if (!multifile.enties && is.data.frame(sheet_entry)) {
          multi.action(paste0(
            "Found ambiguous entry for ID < ", id, " >",
            "\n Matched      : < ", paste0(colnames(sheet)[matches], " > IN:\n", .get.sheet.info.str(sheet)),
            "\n First Match  : ", prev_matches
          ))
          multi.action("Selected the FIRST match !\n")
        } else {
          prev_matches <- paste0("< ", colnames(sheet)[matches], " > IN:\n", .get.sheet.info.str(sheet))
          if (multifile.enties)
            sheet_entry <- append(sheet_entry, list(sheet))
          else
            sheet_entry[[1]] <- sheet
        }
      }
    }

    if (length(sheet_entry) == 0) {
      stop(paste("Did not find entry <", id, ">"))
    }

    # sheet entry is now a list of data.frames (sheets) that hold information about
    # one of more molecules

    # extract time and molecule columns
    mol.infos <- list()
    for (mol in mols) {
      match <- NULL
      for (sheet in sheet_entry) {
        match_tmp <- .id.molecule.cols(sheet, mol$file.name.match,
                                       mol$add.file.matcher, is.fraction = mol$is.fraction,
                                       fixed.unit = mol$fixed.unit, silent = TRUE)
        # found the same entry in multiple sheets
        if (length(match_tmp) > 0 && length(match) > 0)
          stop(paste("Ambiguous sheets found for tag <", mol$file.name.match, ">"))

        # found the same entry in multiple columns
        if (length(match_tmp) > 0 && length(match_tmp$cols) > 1) {
          stop(paste("Ambiguous column found for tag <", mol$file.name.match, ">"))
        }

        if (length(match_tmp) > 0) {
          match <- match_tmp
          match[["sheet"]] <- sheet
        }
      }

      if (length(match) == 0)
        stop(paste("Could not find entry for tag <", mol$file.name.match, ">"))

      #
      #  result <- list(
      #  id = match.tag,
      #  cols = col.nrs,
      #  unit = unit)
      #
      mol.infos <- c(mol.infos, list(match))
    }

    # gather results for each molecule
    pro.results <- list()
    for (j in 1:length(mols)) {
      mol <- mols[[j]]
      mol.i <- mol.infos[[j]]

      time.info <- .id.time.col(mol.i$sheet)
      data <- data.frame(mol.i$sheet[, time.info$col], mol.i$sheet[, mol.i$cols], NA, NA)
      colnames(data) <- c("Time", "Avg", "Min", "Max")

      dup.times <- duplicated(data$Time)
      if (any(dup.times)) {
        message(paste(
          "Deleted duplicated time entries for ID <",
          id, "> and molecule id <", mol$id, ">"
        ))
        data <- data[!dup.times, ]
      }


      entry <- list(
        molecule = mol,
        id = "",
        reference = "",
        group = NA,
        time.unit = time.info$unit,
        value.unit = mol.infos[[j]]$unit,
        data = data,
        type = "individual",
        data.type = "mean",
        origin = "sim"
      )
      class(entry) <- "profile"
      pro.results <- c(pro.results, list(entry))
    }

    plot.infos <- list(
      x.min = master$x.min[i],
      x.max = master$x.max[i],
      y.min = master$y.min[i],
      y.max = master$y.max[i],
      y.min.log = master$y.min.log[i],
      y.max.log = master$y.max.log[i],
      x.offset = master$x.offset[i],
      main = master$main[i]
    )

    sim <- list(
      id = id, obs.ids = obs.ids, profiles = pro.results,
      plot.infos = plot.infos
    )
    class(sim) <- "MatchedProfiles"
    results <- append(results, list(sim))
  }

  return(results)
}
