
.plot.polygon <- function(xvalues, lower, upper, col) {

  x.vals <- c(rev(xvalues), xvalues)
  y.vals <- c(rev(lower), upper)
  graphics::polygon(x.vals, y.vals, border = NA, col = col)
}


plot.profile <- function(profile,
                         poly.color = NA,
                         poly.alpha = 0.5,
                         log = FALSE,
                         ...) {

  add.args <- list(...)
  if (is.na(poly.color)) {
    col <- add.args$col
    if (is.null(col))
      col <- "black"
    if (!is.na(poly.alpha))
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
    if (is.na(profile$data$Min[1])) {
      graphics::points(profile$data$Time, profile$data$Avg, ...)
    } else {
      min.data <- profile$data$Min
      # handle negative mins for log-plots
      if (log) {
        neg.idx <- which(min.data <= 0)
        if (length(neg.idx) > 0) {
          min.data[neg.idx] <- NA
        }
      }

      suppressWarnings(
        gplots::plotCI(x = profile$data$Time, y = profile$data$Avg, gap = 0.0,
            li = min.data, ui = profile$data$Max, add = TRUE, ...))
    }
  } else {
    stop("Not implemented")
  }
}

.legend.args <- function(profile.list) {

  # find the number of unique molecules
  mols <- Map(function(x) if (x$molecule$in.legend) x$molecule else NULL, profile.list)
  mols <- mols[lengths(mols) != 0]
  unique_ids <- unlist(unique(Map(function(x) x$id, mols)))
  mols <- Map(function(x) for (m in mols) if (m$id == x) return(m) , unique_ids)

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
    sim_ref <- Reduce(function(x,y) paste0(x, if (y$molecule$id == tmp.mol$id && y$origin == "sim") y$reference else ""), profile.list, "")
    if (!has_sim)
      mol.ltys[mol_idx] <- NA

    if (nchar(obs_ref) == 0)
      mol.pch[mol_idx] <- NA
    else
      mol.names[mol_idx] <- paste0(mol.names[mol_idx], ", ", obs_ref)

    if (nchar(obs_ref) == 0 && nchar(sim_ref)) {
      mol.names[mol_idx] <- paste0(mol.names[mol_idx], ", ", sim_ref)
    }

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
                         pretty.x.breaks = F,
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

  has_obs <- (length(obs) > 0) && sum(sapply(obs, length)) > 0

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

  if (has_obs)
    all.profiles <- append(obs, matched$profiles)
  else
    all.profiles <- matched$profiles

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
      # observed data with min-max can be plotted event with negative mins
      if (all.profiles[[i]]$origin == "obs") {
        tmp <- all.profiles[[i]]$data[1]
        empty.idx <- which(tmp[-1] <= 0)
        if (length(empty.idx) > 0) {
          all.profiles[[i]]$data <- all.profiles[[i]]$data[-empty.idx,]
        }
      }
      else if (all.profiles[[i]]$data.type == "individual" ) {
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
                      xaxt = if (pretty.x.breaks) "n" else "s",
                      yaxt = if (is.log) "n" else "s",
                      main = main)
  plot.params <- append(plot.params, main.plot.args)
  do.call(graphics::plot, plot.params)
  if (is.log) {
    axis.parameter <- list(pos = "y", range = ylim)
    axis.parameter <- append(axis.parameter, main.plot.args)
    do.call(.minor.tick.log.axis, axis.parameter)
  }

  if (pretty.x.breaks) {
    pb <- 10
    x_unit <- tolower(units::deparse_unit(time.unit))
    if (x_unit == "sec" || x_unit == "min") {

    } else if (x_unit == "h") {
      if (xlim[2] - xlim[1] >= 10)
        pb <- 6
      if (xlim[2] - xlim[1] >= 48)
        pb <- 12
    }

    axis.parameter <- list(prettybase = pb, side = 1, tcl = -0.5, usepar = T, minorn = -1)
    axis.parameter <- append(axis.parameter, main.plot.args)
    do.call(magicaxis::magaxis, axis.parameter)
  }

  # plot profiles
  for (pro in all.profiles) {
    if (pro$origin == "sim")
      lwd = sim.lwd
    else
      lwd = error.lwd

    graphics::plot(pro, col = pro$molecule$color,
                   pch = pro$molecule$pch, lty = pro$molecule$lty, lwd = lwd,
                   log = is.log, ...)
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


.legend_lb <- function(text, col, pch) {
  te <- c()
  co <- c()
  p <- c()
  for (i in 1:length(text)) {
    split <- unlist(strsplit(as.character(text[i]), "\n"))
    l <- length(split)
    te <- c(te, split)
    co <- c(co, col[i], rep(NA_character_, l - 1))
    p <- c(p, pch[i], rep(NA, l - 1))
  }

  return(list(text = te, col = co, pch = p))
}


gof.plot <- function(pred.obs.data, value.lab,
                     target, lwd = 1,
                     group_by = "group",
                     unit.lab = NULL,
                     col = NA,
                     pch = NA,
                     min = NA,
                     max = NA,
                     nice.min = T,
                     nice.max = T,
                     symmetry = F,
                     smart_tagging = F,
                     split_at_group_tag = NA,
                     legend.cex = 1.25,
                     legend.ncol = 1,
                     legend.titles = NULL,
                     show.prop = FALSE,
                     ...) {

  target.pred <- paste0("pred.", target)
  target.obs <- paste0("obs.", target)

  data.pred <- pred.obs.data[target.pred]
  data.obs <- pred.obs.data[target.obs]

  df <- cbind(data.pred, data.obs)
  df <- units::drop_units(df)


  range <- range(df)
  if (!is.na(min))
    range[1] = min

  if (!is.na(max))
    range[2] = max


  if (nice.min) {
    range[1] <- 10^floor(log10(range[1]))
  }

  if (nice.max) {
    range[2] <- 10^ceiling(log10(range[2]))
  }

  if (symmetry) {
    r.log.min <- log10(range[1])
    r.log.max <- log10(range[2])

    log.min <- log10(min(df))
    log.max <- log10(max(df))

    l <- abs(r.log.min - log.min)
    r <- abs(r.log.max - log.max)

    if (l <= 0.25 * r)
      range[1] <- 10^(r.log.min - 1)
    else if (r <= 0.25 * l)
      range[2] <- 10^(r.log.max + 1)
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

  if (show.prop) {
    graphics::curve(1 + 2*(x -1), add = TRUE, from= 1, to = range[2] * 10, lwd=lwd)
    graphics::curve(1/(1 + 2*(1/x -1)), add = TRUE, to= 1, from = range[1] / 10, lwd=lwd)
    graphics::curve(x**2 / (1 + 2*(x -1)), add = TRUE, from= 1, to = range[2] * 10, lwd=lwd)
    graphics::curve(1/((1/x)**2 / (1 + 2*(1/x -1))), add = TRUE, to = 1, from = range[1] / 10, lwd=lwd)
  }

  df <- cbind(df, pred.obs.data[group_by])
  colnames(df) <- c("pred", "obs", "group")
  groups <- unique(df$group)
  n.groups <- length(groups)

  if (length(col) <= 1 && is.na(col))
    col <- seq(from = 1, length.out = n.groups)

  if (length(pch) <= 1 && is.na(pch))
    pch <- seq(from = 15, length.out = n.groups)


  if (smart_tagging) {
    groups <- as.character(groups[order(as.character(groups))])
    tmp_g <- trimws(gsub("\\(.*", "", groups))
    new_col <- col
    new_pch <- pch
    u_g <- unique(tmp_g)
    for (i in 1:length(u_g)) {
      idx <- which(grepl(u_g[i], groups))
      new_col[idx] <- col[i]
      new_pch[idx] <- pch[1:length(idx)]
    }
    col <- new_col
    pch <- new_pch
  }

  i <- 1
  cex <- unlist(list(...)["cex"])
  cex <- if (is.null(cex)) 1.0 else cex
  for (group in groups) {
    subset <- df[df$group == group,]
    graphics::points(subset$obs, subset$pred, pch = pch[i], col = col[i], cex = cex, lwd = cex)
    i <- i + 1
  }

  if (!is.na(split_at_group_tag)) {
    sp_idx <- grep(split_at_group_tag, groups)
    if (length(sp_idx) > 0) {
      col_split <- col[sp_idx]
      pch_split <- pch[sp_idx]
      gr_split <- groups[sp_idx]

      res <- .legend_lb(gr_split, col_split, pch_split)
      col_split <- res$col
      gr_split <- res$text
      pch_split <- res$pch

      if (!is.null(legend.titles[2])) {
        col_split <- c(NA, col_split)
        pch_split <- c(NA, pch_split)
        gr_split <- c(bquote(bold(.(legend.titles[2]))), trimws(gr_split))
      }

      graphics::legend("bottomright", col = col_split, pch = pch_split, legend = as.expression(gr_split),
                       bty = "n",cex = legend.cex, ncol = legend.ncol)

      col <- col[-sp_idx]
      pch <- pch[-sp_idx]
      groups <- groups[-sp_idx]
    }
  }

  res <- .legend_lb(groups, col, pch)
  col <- res$col
  groups <- res$text
  pch <- res$pch

  if (!is.null(legend.titles[1])) {
    col <- c(NA, col)
    pch <- c(NA, pch)
    groups <- c(bquote(bold(.(legend.titles[1]))), trimws(groups))
  }

  graphics::legend("topleft", col = col, pch = pch, legend = as.expression(groups),
                   bty = "n",cex = legend.cex, ncol = legend.ncol)
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
                     split_at_group_tag = NA,
                     legend.ncol = 1,
                     auto.grouping = TRUE,
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

  cex <- unlist(list(...)["cex"])
  cex <- if (is.null(cex)) 1.0 else cex
  if (auto.grouping) {
    groups <- as.character(groups[order(as.character(groups))])
  }

  new_col <- col
  new_pch <- pch
  if (auto.grouping) {
    tmp_g <- trimws(gsub("\\(.*", "", groups))
    u_g <- unique(tmp_g)
    for (i in 1:length(u_g)) {
      idx <- which(grepl(u_g[i], groups))
      new_col[idx] <- col[i]
      new_pch[idx] <- pch[1:length(idx)]
    }
  }
  col <- new_col
  pch <- new_pch

  i <- 1
  for (group in groups) {
    subset <- data %>% dplyr::filter(.[[group_by]] == group)
    graphics::points(subset$Obs, subset$Pred, pch = pch[i], col = col[i], cex = cex, lwd = cex)
    i <- i + 1
  }

  if (!is.na(split_at_group_tag)) {
    sp_idx <- grep(split_at_group_tag, groups)
    if (length(sp_idx) > 0) {
      col_split <- col[sp_idx]
      pch_split <- pch[sp_idx]
      gr_split <- groups[sp_idx]

      res <- .legend_lb(gr_split, col_split, pch_split)
      col_split <- res$col
      gr_split <- res$text
      pch_split <- res$pch

      graphics::legend("bottomright", col = col_split, pch = pch_split, legend = trimws(gr_split),
                       bty = "n",cex = legend.cex, ncol = legend.ncol)

      col <- col[-sp_idx]
      pch <- pch[-sp_idx]
      groups <- groups[-sp_idx]
    }
  }

  res <- .legend_lb(groups, col, pch)
  col <- res$col
  groups <- res$text
  pch <- res$pch

  graphics::legend("topleft", col = col, pch = pch, legend = trimws(groups),
                   bty = "n",cex = legend.cex, ncol = legend.ncol)

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


