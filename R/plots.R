
.plot.polygon <- function(xvalues, lower, upper, col) {

  x.vals <- c(rev(xvalues), xvalues)
  y.vals <- c(rev(lower), upper)
  graphics::polygon(x.vals, y.vals, border = NA, col = col)
}


plot.profile <- function(profile,
                         poly.color = NA,
                         poly.alpha = 0.5, ...) {

  add.args <- list(...)
  if (is.na(poly.color)) {
    col <- add.args$col
    if (is.null(col))
      col <- "black"
    poly.color <- add.alpha(col, poly.alpha)
  }

  if (!is.profile(profile))
    stop("Input must be of class profile")

  if (profile$origin == "sim" && profile$data.type == "mean") {
    .plot.polygon(profile$data$Time, profile$data$Min, profile$data$Max, col = poly.color)
    graphics::lines(profile$data$Time, profile$data$Avg, ...)

  } else if (profile$origin == "sim" && profile$data.type == "individual") {
      data <- profile$data[-1]
      for (i in 1:ncol(data)) {
        graphics::lines(profile$data$Time, data[,i], ...)
    }

  } else if (profile$origin == "obs" && profile$data.type == "mean") {
    if (is.na(profile$data$Min[1]))
      graphics::points(profile$data$Time, profile$data$Avg, ...)
    else
      suppressWarnings(
        gplots::plotCI(x = profile$data$Time, y = profile$data$Avg, gap = 0.0,
            li = profile$data$Min, ui = profile$data$Max, add = TRUE, ...))

  } else {
    stop("Not implemented")
  }
}

.legend.args <- function(profile.list) {

  # find the number of unique molecules
  mols <- Map(function(x) x$molecule, profile.list)
  unique_ids <- unlist(unique(Map(function(x) x$id, mols)))
  mols <- Map(function(x) for(m in mols) if (m$id == x) return(m) , unique_ids)

  mol.names <- Reduce(function(x,y) append(x, y$display.name), mols, c())
  mol.cols <- Reduce(function(x,y) {append(x, y$color)}, mols, c())
  mol.ltys <- Reduce(function(x,y) {append(x, y$lty)}, mols, c())
  mol.pch <- Reduce(function(x,y) {append(x, y$pch)}, mols, c())

  # points only for observed data
  # lines only for simulated data
  for (mol_idx in 1:length(mols)) {
    tmp.mol <- mols[[mol_idx]]
    has_sim <- Reduce(function(x,y) x || (y$molecule$id == tmp.mol$id && y$origin == "sim"), profile.list, FALSE)
    obs_ref <- Reduce(function(x,y) paste0(x, if (y$molecule$id == tmp.mol$id && y$origin == "obs") y$reference else "") , profile.list, "")
    if (!has_sim)
      mol.ltys[mol_idx] <- NA

    if (nchar(obs_ref) == 0)
      mol.pch[mol_idx] <- NA
    else
      mol.names[mol_idx] <- paste0(mol.names[mol_idx], ", ", obs_ref)
  }

  # sort by name
  df <- data.frame(mol.names, mol.cols, mol.ltys, mol.pch, stringsAsFactors = F)
  df <- df[order(mol.names),]

  return(list(legend = df$mol.names, col = df$mol.cols, lty = df$mol.ltys, pch = df$mol.pch))
}





plot.matched <- function(matched, obs.data,
                         time.unit = NA,
                         value.unit = NA,
                         ylab = NA,
                         xlab = "Time",
                         xlim = NA,
                         ylim = NA,
                         avg.fn = base::mean,
                         min.var.fn = std.dev.min,
                         max.var.fn = std.dev.max,
                         rm.zero.neg.rows = F,
                         ymax.rel.add = 0.125,
                         show.legend = T,
                         show.main = T,
                         sim.lwd = 2.8,
                         error.lwd = 1.8,
                         legend.plot.args = list(x = "topright", cex = 1.25, lwd = 3, bty = "n"),
                         main.plot.args = list(bty = 'l', las = 1, cex.axis = 1.5, cex.lab = 1.7, cex = 1.5),
                         ...) {

  if (!is.matched.profiles(matched))
    stop("Input must be of class matchedProfiles")

  # gather observed data
  obs <- .gather.ids(obs.data, matched$obs.ids)
  if (length(obs) == 0)
    message(paste("Attention: Matched profiles with id < ", matched$id ,"> do not have any observed data"))

  if (!is.na(matched$obs.ids) && length(obs) != length(matched$obs.ids))
    stop(paste("Attention: Matched profiles with id < ", matched$id ,"> have missing observed data"))

  has_obs <- length(obs) > 0

  # average profile data
  for (i in 1:length(matched$profiles)) {
      if (rm.zero.neg.rows) {
        if (matched$profiles[[i]]$data.type == "individual") {
          tmp <- matched$profiles[[i]]$data
          matched$profiles[[i]]$data <- tmp[rowSums(tmp[-1] <= 0 ) == 0, ]
        } else {
          tmp <- matched$profiles[[i]]$data
          row_sub = apply(tmp, 1, function(row) all(row > 0, na.rm = T))
          matched$profiles[[i]]$data <- tmp[row_sub,]
        }
      }

      if (matched$profiles[[i]]$data.type == "individual") {
        matched$profiles[[i]] <- average.pop.profile(matched$profiles[[i]],
                                                     avg.fn = avg.fn,
                                                     min.var.fn = min.var.fn,
                                                     max.var.fn = max.var.fn)
      }
  }

  # match units
  if (is.na(time.unit)) {
    if (has_obs)
      time.unit <- obs[[1]]$time.unit
    else
      time.unit <- matched$profiles[[1]]$time.unit
  }

  if (is.na(value.unit)) {
    if (has_obs)
      value.unit <- obs[[1]]$value.unit
    else
      value.unit <- matched$profiles[[1]]$value.unit
  }

  all.profiles <- append(obs, matched$profiles)
  for (i in 1:length(all.profiles)) {
    all.profiles[[i]] <- convert.profile(all.profiles[[i]], value.unit = value.unit, time.unit = time.unit)
  }

  # new x-limits -> trim the data
  if (!is.blank(matched$plot.infos$x.max) && !is.blank(matched$plot.infos$x.min)) {
    x.min <- matched$plot.infos$x.min
    units(x.min) <- time.unit
    x.min <- units::drop_units(x.min)
    x.max <- matched$plot.infos$x.max
    units(x.max) <- time.unit
    x.max <- units::drop_units(x.max)

    for (i in 1:length(all.profiles)) {
      all.profiles[[i]] <- trim.time(all.profiles[[i]], from = x.min, to = x.max)
    }
  }

  is.log <-  ("log" %in% names(main.plot.args) &&
                grepl("y", main.plot.args[["log"]]))

  if (is.log) {
    if (!is.blank(matched$plot.infos$y.min.log))
      matched$plot.infos$y.min <- matched$plot.infos$y.min.log

    if (!is.blank(matched$plot.infos$y.max.log))
      matched$plot.infos$y.max <- matched$plot.infos$y.max.log


    for (i in 1:length(all.profiles)) {
      if (all.profiles[[i]]$data.type == "individual") {
        tmp <- all.profiles[[i]]$data
        all.profiles[[i]]$data <- tmp[rowSums(tmp[-1] <= 0 ) == 0, ]
      } else {
        tmp <- all.profiles[[i]]$data
        row_sub = apply(tmp, 1, function(row) all(row > 0, na.rm = T))
        all.profiles[[i]]$data <- tmp[row_sub,]
      }
    }
  }

  if (!is.blank(matched$plot.infos$x.offset)) {
    x_offset <- matched$plot.infos$x.offset
    units(x_offset) <- time.unit
    x_offset <- units::drop_units(x_offset)

    for (i in 1:length(all.profiles)) {
      all.profiles[[i]]$data$Time <- all.profiles[[i]]$data$Time - x_offset
    }
  }

  y.min <- NA
  y.max <- NA
  if (!is.blank(matched$plot.infos$y.max) && !is.blank(matched$plot.infos$y.min)) {
    y.min <- matched$plot.infos$y.min
    units(y.min) <- value.unit
    y.min <- units::drop_units(y.min)
    y.max <- matched$plot.infos$y.max
    units(y.max) <- value.unit
    y.max <- units::drop_units(y.max)
  }

  if (!is.na(y.min) && !is.na(y.max))
    ylim <- c(y.min, y.max)

  # ranges
  if (length(xlim) <= 1 && is.na(xlim)) {
    xlim <- c(Inf, -Inf)
    for (pro in all.profiles) {
      x.range <- range(pro, range.type = "time")
      xlim[1] <- min(xlim[1], x.range[1])
      xlim[2] <- max(xlim[2], x.range[2])
    }
  }

  if (length(ylim) <= 1 && is.na(ylim)) {
    ylim <- c(Inf, -Inf)
    for (pro in all.profiles) {
      y.range <- range(pro, range.type = "data")
      ylim[1] <- min(ylim[1], y.range[1])
      ylim[2] <- max(ylim[2], y.range[2])
    }
  }

  if (is.log && ylim[1] <= 0) {
    ylim[1] <- 0 + ylim[2]/1000
  }

  # factor
  ylim[2] <- ylim[2] + abs(ylim[1] - ylim[2]) * ymax.rel.add

  # labs
  units::units_options(parse = FALSE)
  xlab <- units::make_unit_label(xlab, time.unit)
  if (is.na(ylab)) {
    ylab <- all.profiles[[1]]$molecule$ylab
  }
  ylab <- units::make_unit_label(ylab, value.unit)
  units::units_options(parse = TRUE)

  # check for main
  main = NA
  if (show.main)
    main <- matched$plot.infos$main

  # HACK !!
  if (rm.zero.neg.rows && ylim[1] < 0) {
    ylim[1] = 1E-4
  }

  plot.params <- list(1, type = "n", xlim = xlim, ylim = ylim,
                      xlab = xlab,
                      ylab = ylab,
                      yaxt = if (is.log) "n" else "s",
                      main = main)
  plot.params <- append(plot.params, main.plot.args)
  do.call(graphics::plot, plot.params)
  if (is.log) {
    axis.parameter <- list(pos = "y", range = ylim)
    axis.parameter <- append(axis.parameter, main.plot.args)
    do.call(.minor.tick.log.axis, axis.parameter)
  }

  # plot profiles
  for (pro in all.profiles) {
    if (pro$origin == "sim")
      lwd = sim.lwd
    else
      lwd = error.lwd

    graphics::plot(pro, col = pro$molecule$color, pch = pro$molecule$pch, lty = pro$molecule$lty, lwd = lwd, ...)
  }


  if (show.legend) {
    leg.data <- .legend.args(all.profiles)
    do.call(graphics::legend, append(legend.plot.args, leg.data))
  }
}


.minor.tick.log.axis <- function(pos = c("x", "y"), range, ...) {

  pos <- match.arg(pos)

  y1  <- floor(log10(range))
  pow <- seq(y1[1], y1[2] + 1)
  pow.labels <- sapply(pow,function(i)
    as.expression(bquote(10^ .(i)))
  )

  ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))

  pos.id = 2
  if (pos == "x")
    pos.id = 1

  graphics::axis(pos.id, 10^pow, las = 1, labels = pow.labels, ...)
  graphics::axis(pos.id, ticksat, labels = NA, tcl = -0.25, lwd = 0, lwd.ticks = 1, ...)
}


gof.plot <- function(pred.obs.data, value.lab,
                     target, lwd = 1,
                     group_by = "group",
                     unit.lab = NULL,
                     col = NA,
                     pch = NA,
                     nice.min = T,
                     nice.max = T,
                     legend.cex = 1.25,
                     ...) {

  target.pred <- paste0("pred.", target)
  target.obs <- paste0("obs.", target)

  data.pred <- pred.obs.data[target.pred]
  data.obs <- pred.obs.data[target.obs]

  df <- cbind(data.pred, data.obs)
  df <- units::drop_units(df)


  range <- range(df)
  if (nice.min) {
    range[1] <- 10^floor(log10(range[1]))
  }

  if (nice.max) {
    range[2] <- 10^ceiling(log10(range[2]))
  }

  unit <- units(data.pred[,1])

  lab_unit <- if (is.null(unit.lab)) units::make_unit_label("", unit, parse = T) else unit.lab
  ylab <- substitute(a ~ b ~ c, lapply(list(a = "Predicted",
                                            b = value.lab,
                                            c = lab_unit), "[[", 1))

  xlab <- substitute(a ~ b ~ c, lapply(list(a = "Observed",
                                            b = value.lab,
                                            c = lab_unit), "[[", 1))
  # base plot
  graphics::plot(df, log = "xy", type = "n",
                 ylab = ylab,
                 xlab = xlab,
                 xlim = range,
                 ylim = range,
                 yaxt = "n",
                 xaxt = "n",
                 ...)
  .minor.tick.log.axis("x", range, ...)
  .minor.tick.log.axis("y", range, ...)

  # ident line plus 2- and 1.25-fold range
  graphics::abline(0, 1, lwd = lwd)
  graphics::abline(0 , 2, lty = 2, untf = T, lwd = lwd)
  graphics::abline(0 , 0.5, lty = 2, untf = T, lwd = lwd)
  graphics::abline(0, 1.25, lty = 3, untf = T, lwd = lwd)
  graphics::abline(0, 1/1.25, lty = 3, untf = T, lwd = lwd)

  df <- cbind(df, pred.obs.data[group_by])
  colnames(df) <- c("pred", "obs", "group")
  groups <- unique(df$group)
  n.groups <- length(groups)

  if (length(col) <= 1 && is.na(col))
    col <- seq(from = 1, length.out = n.groups)

  if (length(pch) <= 1 && is.na(pch))
    pch <- seq(from = 15, length.out = n.groups)

  i <- 1
  cex <- unlist(list(...)["cex"])
  cex <- if (is.null(cex)) 1.0 else cex
  for (group in groups) {
    subset <- df[df$group == group,]
    graphics::points(subset$obs, subset$pred, pch = pch[i], col = col[i], cex = cex)
    i <- i + 1
  }

  graphics::legend("topleft", col = col, pch = pch, legend = groups,
                   bty = "n",cex = legend.cex)
}

pred_obs_plot <- function(pred.obs.data,
                     value.lab = "Concentration",
                     lwd = 1,
                     group_by = "ref",
                     match_by = NA,
                     match = NA,
                     col = NA,
                     pch = NA,
                     nice.min = T,
                     nice.max = T,
                     legend.cex = 1.25,
                     show.two.fold = F,
                     show.1.25.fold = F,
                     ...) {

  data <- pred.obs.data$data
  if (!is.na(match[1]) && !is.na(match_by[1])) {
    data <- data %>% dplyr::filter(.[[match_by]] == match)
  }

  if (nrow(data) == 0) {
    stop("pred.obs.data has not data attached or filter did not match")
  }

  data <- dplyr::select(data, "Pred", "Obs", group_by)
  range <- range(data[1:2])
  if (nice.min) {
    range[1] <- 10^floor(log10(range[1]))
  }

  if (nice.max) {
    range[2] <- 10^ceiling(log10(range[2]))
  }

  unit <- pred.obs.data$meta$value.unit

  lab_unit <- units::make_unit_label("", unit, parse = T)
  ylab <- substitute(a ~ b ~ c, lapply(list(a = "Predicted",
                                            b = value.lab,
                                            c = lab_unit), "[[", 1))

  xlab <- substitute(a ~ b ~ c, lapply(list(a = "Observed",
                                            b = value.lab,
                                            c = lab_unit), "[[", 1))
  # base plot
  graphics::plot(1, log = "xy", type = "n",
                 ylab = ylab,
                 xlab = xlab,
                 xlim = range,
                 ylim = range,
                 yaxt = "n",
                 xaxt = "n",
                 ...)
  .minor.tick.log.axis("x", range, ...)
  .minor.tick.log.axis("y", range, ...)

  # ident line plus 2- and 1.25-fold range (optional)
  graphics::abline(0, 1, lwd = lwd)
  if (show.two.fold) {
    graphics::abline(0 , 2, lty = 2, untf = T, lwd = lwd)
    graphics::abline(0 , 0.5, lty = 2, untf = T, lwd = lwd)
  }

  if (show.1.25.fold) {
    graphics::abline(0, 1.25, lty = 3, untf = T, lwd = lwd)
    graphics::abline(0, 1/1.25, lty = 3, untf = T, lwd = lwd)
  }

  groups <- dplyr::pull(unique(data[3]))
  n.groups <- length(groups)

  if (length(col) <= 1 && is.na(col))
    col <- seq(from = 1, length.out = n.groups)

  if (length(pch) <= 1 && is.na(pch))
    pch <- seq(from = 15, length.out = n.groups)

  i <- 1
  cex <- unlist(list(...)["cex"])
  cex <- if (is.null(cex)) 1.0 else cex
  for (group in groups) {
    subset <- data %>% dplyr::filter(.[[group_by]] == group)
    graphics::points(subset$Obs, subset$Pred, pch = pch[i], col = col[i], cex = cex)
    i <- i + 1
  }

  graphics::legend("topleft", col = col, pch = pch, legend = groups,
                   bty = "n",cex = legend.cex)
}


axis.labels <- function(profile,
                        x.prefix = "Time",
                        y.prefix = "Concentration") {

  if (!is.profile(profile))
    stop("Input must be of class profile")

  x <- units::make_unit_label(x.prefix, profile$time.unit, parse = T)
  y <- units::make_unit_label(y.prefix, profile$value.unit, parse = T)

  return(c(x,y))
}


