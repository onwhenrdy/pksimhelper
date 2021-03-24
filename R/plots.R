
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

      dots <- list(...)
      dots$lty = 1
      args = list(x = profile$data$Time, y = profile$data$Avg, gap = 0.0,
                  li = min.data, ui = profile$data$Max, add = TRUE)
      args <- c(args, dots)
      suppressWarnings(do.call(gplots::plotCI, args))
    }
  } else {
    stop("Not implemented")
  }
}

.legend.args <- function(profile.list, format = "{mol}, {ref}") {

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
    obs_ref <- Reduce(function(x,y) paste0(x, if (y$molecule$id == tmp.mol$id && y$origin == "obs" && nchar(x) == 0) y$reference else "") , profile.list, "")
    sim_ref <- Reduce(function(x,y) paste0(x, if (y$molecule$id == tmp.mol$id && y$origin == "sim") y$reference else ""), profile.list, "")
    
    N <- Reduce(function(x,y) paste0(x, if (y$molecule$id == tmp.mol$id && y$origin == "obs") y$N else ""), profile.list, "")
    
    if (!has_sim)
      mol.ltys[mol_idx] <- NA

    if (nchar(obs_ref) == 0)
      mol.pch[mol_idx] <- NA
    else {
      
      mol <- mol.names[mol_idx]
      ref <- obs_ref
      mol.names[mol_idx] <- glue::glue(format)
    }
    
    if (nchar(obs_ref) == 0 && nchar(sim_ref)) {
      
      mol <- mol.names[mol_idx]
      ref <- sim_ref
      mol.names[mol_idx] <- glue::glue(format)
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
                         add.legend.text = NA,
                         legend.plot.args = list(x = "topright", cex = 1.25, lwd = 3, bty = "n"),
                         main.plot.args = list(bty = 'l', las = 1, cex.axis = 1.5, cex.lab = 1.7, cex = 1.5),
                         legend_format = "{mol}, {ref}",
                         ...) {

  if (!is.matched.profiles(matched))
    stop("Input must be of class matchedProfiles")

  # gather observed data
  obs <- .gather.ids(obs.data, matched$obs.ids)
  if (length(obs) == 0)
    message(paste("Attention: Matched profiles with id < ", matched$id ,"> do not have any observed data"))

  if (!is.na(matched$obs.ids) && length(obs) != length(matched$obs.ids))
    stop(paste("Attention: Matched profiles with id < ", matched$id ,"> have missing observed data"))

  obs <- purrr::compact(obs)
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

  
  # LOG
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

  
  # x-offset
  offset <- function(x_offsets, all_profiles) {
    if(length(x_offsets) < length(all_profiles))
      x_offsets <- rep(x_offsets, length(all_profiles))
    
    for (i in 1:length(all_profiles)) {
      x_offset <- x_offsets[i]
      units(x_offset) <- time.unit
      x_offset <- units::drop_units(x_offset)
      
      all_profiles[[i]]$data$Time <- all_profiles[[i]]$data$Time - x_offset
    }
    
    return(all_profiles)
  }
  
  if (!is.blank(matched$plot.infos$x.offset))
  {
    # simple version
    if (!is.list(matched$plot.infos$x.offset))
      all.profiles <- offset(matched$plot.infos$x.offset, all.profiles)
    else {
      
      obs_n <- if (!has_obs) 0 else length(obs)
      if (has_obs) {
        all.profiles[1:obs_n] <- offset(matched$plot.infos$x.offset$obs, all.profiles[1:obs_n])
      }
      
      all.profiles[(obs_n+1):length(all.profiles)] <- offset(matched$plot.infos$x.offset$profiles, 
                                                       all.profiles[(obs_n+1):length(all.profiles)])
    }
  }

  # new x-limits -> trim the data
  if (!is.blank(matched$plot.infos$x.max) || !is.blank(matched$plot.infos$x.min)) {
    
    x.min <- NA
    if (!is.blank(matched$plot.infos$x.min)) {
      x.min <- matched$plot.infos$x.min
      units(x.min) <- time.unit
      x.min <- units::drop_units(x.min)
    }
    
    x.max <- NA
    if (!is.blank(matched$plot.infos$x.max)) {
      x.max <- matched$plot.infos$x.max
      units(x.max) <- time.unit
      x.max <- units::drop_units(x.max)
    }
    
    for (i in 1:length(all.profiles)) {
      all.profiles[[i]] <- trim.time(all.profiles[[i]], from = x.min, to = x.max)
    }
    
  }
  
  # y-limits
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
    leg.data <- .legend.args(all.profiles, format = legend_format)
    if (!is.na(add.legend.text)) {
      leg.data$legend <- c(leg.data$legend, add.legend.text)
      leg.data$col <- c(leg.data$col, NA)
      leg.data$lty <- c(leg.data$lty, NA)
      leg.data$pch <- c(leg.data$pch, NA)
    }
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

  suppressWarnings(graphics::axis(pos.id, 10^pow, las = 1, labels = pow.labels, ...))
  suppressWarnings(graphics::axis(pos.id, ticksat, labels = NA, tcl = -0.25, 
                                  lwd = 0, lwd.ticks = 1, ...))
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


# Deprecated
gof.plot <- function(pred.obs.data, value.lab,
                     target, lwd = 1,
                     group_by = "group",
                     subgroup_by = NA,
                     unit.lab = NULL,
                     col = NA,
                     pch = NA,
                     min = NA,
                     max = NA,
                     nice.min = T,
                     nice.max = T,
                     symmetry = F,
                     smart_tagging = F,
                     smart_tagging_pattern = "\\(.*",
                     split_at_group_tag = NA,
                     legend.cex = 1.25,
                     legend.ncol = 1,
                     legend.titles = NULL,
                     show.prop = FALSE,
                     ...) {

  .Deprecated("plot_gof_pk")
  
  
  if (smart_tagging && !is.na(subgroup_by)) {
    stop("smart_tagging and subgroup_by options are not compatible")
  }


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
  if (!is.na(subgroup_by)) {
    df <- cbind(df, pred.obs.data[subgroup_by])
    colnames(df) <- c("pred", "obs", "group", "subgroup")
  }

  groups <- unique(df$group)
  n.groups <- length(groups)

  if (length(col) <= 1 && is.na(col))
    col <- seq(from = 1, length.out = n.groups)

  if (length(pch) <= 1 && is.na(pch))
    pch <- seq(from = 15, length.out = n.groups)

  if (smart_tagging) {
    groups <- as.character(groups[order(as.character(groups))])
    tmp_g <- trimws(gsub(smart_tagging_pattern, "", groups))
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

  cex <- unlist(list(...)["cex"])
  cex <- if (is.null(cex)) 1.0 else cex

  new_col <- col
  new_pch <- pch
  ##########################################################################
  tmp_groups <- groups
  if (!is.na(subgroup_by)) {
    col <- c()
    pch <- c()
    subgroups <- as.character(unique(df[["subgroup"]]))
    groups <- c()
    pch_map <- as.list(new_pch[1:length(subgroups)])
    names(pch_map) <- subgroups

    message("Point subgroup mapping:")
    message(paste(paste(" ", names(pch_map)),
                  pch_map,
                  sep = " = ", collapse = "\n"))
  }

  i <- 1
  for (group in tmp_groups) {
    subset <- df[df$group == group,]

    if (!is.na(subgroup_by)) {
      for (sub in subgroups) {

        subset_2 <- subset %>% dplyr::filter(.[["subgroup"]] == sub)
        if (nrow(subset_2) > 0) {
          pch <- c(pch, pch_map[[sub]])
          col <- c(col, new_col[i])
          groups <- c(groups, as.character(subset_2[["group"]][1]))
          graphics::points(subset_2$obs,
                           subset_2$pred,
                           pch = pch_map[[sub]],
                           col = new_col[i],
                           cex = cex,
                           lwd = cex)
        }
      }

      i <- i + 1
    } else  {
      graphics::points(subset$obs,
                       subset$pred,
                       pch = pch[i],
                       col = col[i],
                       cex = cex,
                       lwd = cex)
      i <- i + 1
    }
  }

  if (!is.na(split_at_group_tag)) {
    if (is.numeric(split_at_group_tag)) {
      sp_idx <- split_at_group_tag:length(groups)
    } else {
      sp_idx <- grep(split_at_group_tag, groups)
    }
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

      graphics::legend("bottomright",
                       col = col_split,
                       pch = pch_split,
                       legend = as.expression(gr_split),
                       bty = "n",
                       cex = legend.cex,
                       ncol = legend.ncol)

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
                     subgroup_by = NA,
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

  if (auto.grouping && !is.na(subgroup_by)) {
    stop("auto.grouping and subgroup_by options are not compatible")
  }

  data <- pred.obs.data$data
  if (!is.na(match[1]) && !is.na(match_by[1])) {
    data <- data %>% dplyr::filter(.[[match_by]] == match)
  }

  if (nrow(data) == 0) {
    stop("pred.obs.data has not data attached or filter did not match")
  }

  if (!is.na(subgroup_by))
    data <- dplyr::select(data, "Pred", "Obs", group_by, subgroup_by)
  else
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

  col_df <- col
  if (is.data.frame(col_df) || (length(col) <= 1 && is.na(col)))
    col <- .get_distinct_colors(n.groups)

  if (length(pch) <= 1 && is.na(pch))
    pch <- seq(from = 15, length.out = n.groups)

  cex <- unlist(list(...)["cex"])
  cex <- if (is.null(cex)) 1.0 else cex
  if (auto.grouping) {
    groups <- as.character(groups[order(as.character(groups))])
  }

  if (is.data.frame(col_df)) {
    col <- .extract_color_map_col(col_df, groups)
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



  ##########################################################################
  tmp_groups <- groups
  if (!is.na(subgroup_by)) {
    col <- c()
    pch <- c()
    subgroups <- unique(data[[subgroup_by]])
    groups <- c()
    pch_map <- as.list(new_pch)
    names(pch_map) <- subgroups

    message("* Point subgroup mapping:")
    message(paste(paste(" -", names(pch_map)),
                  pch_map,
                  sep = " => ", collapse = "\n"))
  }

  i <- 1
  for (group in tmp_groups) {

    #subset <- data %>% dplyr::filter(.[[group_by]] == group)
    subset <- data[data[group_by] == group,]
    # this is special handling !!
    if (!is.na(subgroup_by)) {
      for (sub in subgroups) {
        subset_2 <- subset %>% dplyr::filter(.[[subgroup_by]] == sub)
        if (nrow(subset_2) > 0) {
          pch <- c(pch, pch_map[[sub]])
          col <- c(col, new_col[i])
          groups <- c(groups, as.character(subset_2[[group_by]][1]))
          graphics::points(subset_2$Obs,
                           subset_2$Pred,
                           pch = pch_map[[sub]],
                           col = new_col[i],
                           cex = cex, lwd = cex)
        }
      }
      i <- i + 1

    } else  {
      graphics::points(subset$Obs, subset$Pred,
                       pch = pch[i],
                       col = col[i],
                       cex = cex, lwd = cex)
      i <- i + 1
    }
  }

  # Here: Groups, pch and cols must be set !!!
  ####################################################################################

  if (!is.na(split_at_group_tag)) {
    if (is.numeric(split_at_group_tag)) {
      sp_idx <- split_at_group_tag:length(groups)
    } else {
      sp_idx <- grep(split_at_group_tag, groups)
    }
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


# TODO: md_assist
plot_profiles <- function(profiles, observed_data, md_assist = NULL,
                          plot_folder = NULL,
                          format = c("pdf", "tiff", "png", "jpeg"),
                          format_opts = NULL,
                          file_prefix = NULL,
                          counter_prefix = FALSE,
                          x_unit = NA,
                          y_unit = NA,
                          xlab = "Time", 
                          pretty_x_breaks = TRUE,
                          plot_log = TRUE,
                          plot_lin = TRUE,
                          pop_fn = list(avg = geo.mean, min = geo.dev.min, max = geo.dev.max),
                          geom = list(cex = 1.6, lwd = 3, err_lwd = 2.5, poly_alpha = 0.2),
                          par_fn = NULL,
                          plot_args = NULL,
                          legend_args = NULL,
                          legend_format = "{mol}, {ref}",
                          add_legend_text = NULL, # "N/n" or "text/text vector" or function (profile, counter, obs_data)
                          panel_first = NULL,
                          sanatize_id = TRUE,
                          silent = FALSE,
                          parallel = FALSE,
                          parallel.cores = "auto",
                          ...) {
  
  # validate
  format <- match.arg(format)
  
  # prepare functions
  msg_fn <- if (silent) function(x, nl = TRUE) {} else function(x, nl = TRUE) {message(x, appendLF = nl)}
  
  plot_fn <- function(file_name) {
      do.call(format, append(list(file_name), format_opts))
  }
  
  id_fn <- function(profile, counter, is_fraction, linlog = c("linear", "log")) {
    
    # for user functions
    if (is.function(file_prefix)) {
      
      result <- file_prefix(profile, counter, is_fraction, linlog)
      if(!is.character(result))
        stop("file_prefix function does not return a string")  
      return(result)
    }
    
    id <- profile$id
    if (sanatize_id)
      id <- .sanatize_filename(id)
    
    
    # standard is (prefix)_ID_(frac)_(1)[LOG/LIN]
    result <- if (!is.null(file_prefix)) paste0(file_prefix, "_", id) else id
    
    if (is_fraction)
      result <- paste0(result, "_frac")
    
    if (counter_prefix)
      result <- paste0(result, "_", counter)
    
    result <- paste0(result, "_", toupper(linlog))
    
    return(result)
  }
  
  units_fn <- function(profile, counter, is_fraction) {
    
    unit_time <- NA
    unit_value <- NA
    
    if (is.function(x_unit))
      unit_time <- x_unit(profile, counter, is_fraction)
    else
      unit_time <- x_unit
    
    if (is.function(y_unit))
      unit_value <- y_unit(profile, counter, is_fraction)
    else {
      unit_value <- y_unit
      if (!is.na(unit_value) && is_fraction)
        unit_value <- as_units("%")
    }
    
    result <- list(time = unit_time,
                   value = unit_value)
    
    return(result)
  }
  
  par_plot_fn <- function() {
    
    if (!is.null(par_fn))
      par_fn()
    else
      par(mar = par()$mar + c(0, 2, 0, 0), mgp = par()$mgp + c(0.2, -0.0, 0))
  }
  
  panel_fn <- function(profile, counter, is_fraction, linlog) {
    if (!is.null(panel_first))
      panel_first(profile, counter, is_fraction, linlog)
  }
  
  
  plot_lin_fn <- function(profile, counter, is_fraction) {
    if (is.function(plot_lin))
      return(plot_lin(profile, counter, is_fraction))
    
    return(plot_lin)
  }
  
  plot_log_fn <- function(profile, counter, is_fraction) {
    if (is.function(plot_log))
      return(plot_log(profile, counter, is_fraction))
    
    return(plot_log && !is_fraction)
  }
  
  add_legend_fn <- function(profile, counter, obs_data) {
    
    if(is.character(add_legend_text)) {
      
      if (add_legend_text == "n" || add_legend_text == "N") {
        
        # number of subjects
        tmp <- paste(add_legend_text, "=")
        
        # gather observed data
        obs <- .gather.ids(obs_data, profile$obs.ids)
        obs <- purrr::compact(obs)
        if (length(obs) == 0)
          return(NA)
        
        if (!is.na(profile$obs.ids) && length(obs) != length(profile$obs.ids))
          stop(paste("Matched profiles with id < ", profile$id ,"> have missing observed data"))
        

        n_sub <- sapply(obs, function(x) {x$N})
        n_sub <- unique(n_sub)
        if (length(n_sub) > 1) {
          msg_fn(paste("\n  * Found different N for observed data. First non-NA was used."))
          
          n_sub <- dplyr::first(na.omit(n_sub))
        }
        
        return(paste(tmp, n_sub))
        
      } else {
        # single string or vector of strings
        if (length(add_legend_text) > 1)
          return(add_legend_text[counter])
        else
          return(add_legend_text)
      }
      
      
    } else if (is.function(add_legend_text)) {
      
      # custom function
      # gather observed data
      obs <- .gather.ids(obs_data, profile$obs.ids)
      obs <- purrr::compact(obs)
      return(add_legend_text(profile, counter, obs))
    }
    
    return(NA)
  }
  
  
  ###############################################################################
  if (is.null(names(pop_fn)))
    names(pop_fn) <- c("avg", "min", "max")
  avg.fn = pop_fn$avg
  min.var.fn = pop_fn$min
  max.var.fn = pop_fn$max
  
  
  main_plot_args_def <- list(bty = 'l', las = 0, cex.main = 1.6,
                             cex.axis = 1.5, cex.lab = 1.5, 
                             cex = 1.5)
  main_plot_args <- .ls_override(main_plot_args_def, plot_args)
  
  legend_plot_args_def <- list(x = "topright", cex = 1.4, 
                               lwd = 3, bty = "n", 
                               y.intersp = 1.2)
  legend_plot_args <- .ls_override(legend_plot_args_def, legend_args)
  
  cex <- geom$cex
  poly.alpha <- geom$poly_alpha
  sim.lwd <- geom$lwd
  error.lwd <- geom$err_lwd
  
  # plotting loop
  msg_fn(paste("Plotting to folder: <", plot_folder, ">"))
  
  if (parallel) {
    if (is.character(parallel.cores)) {
      if (tolower(parallel.cores) == "auto")
        parallel.cores <- min(parallel::detectCores(logical = FALSE), 
                              max(1, floor(length(profiles) / 2)))
      else
        stop("parallel.cores must be auto or a positive number", call. = FALSE)
    } else {
      parallel.cores <- as.numeric(parallel.cores)
      if (parallel.cores <= 0)
        stop("parallel.cores must be a positive number")
    }
    
    msg_fn(paste("Using <", parallel.cores, "> cores for parallel computing"))
    doParallel::registerDoParallel(cores = parallel.cores)
    on.exit(doParallel::stopImplicitCluster())
  }
  
  `%doit%` <- if(parallel) foreach::`%dopar%` else foreach::`%do%`
  foreach::foreach(i = 1:length(profiles)) %doit% {
    profile <- profiles[[i]]
    counter <- i
    msg_fn(paste("* Processing <", profile$id, ">"), nl = FALSE)
    
    is_fraction <- has_fraction(profile)
    units <- units_fn(profile, counter, is_fraction)
    
    add_led <- add_legend_fn(profile, counter, observed_data)
    
    # linear plots
    p_lin <- plot_lin_fn(profile, counter, is_fraction)
    msg_fn(paste("\tLin:", p_lin), nl = FALSE)
    if (p_lin) {
      
      file_name <- paste0(id_fn(profile, counter, is_fraction, "linear"), ".", format)
      full_name <- file.path(plot_folder, file_name)
      
      plot_fn(full_name)
      par_plot_fn()
      
      plot.matched(profile, observed_data,
                   avg.fn = avg.fn, min.var.fn = min.var.fn, max.var.fn = max.var.fn,
                   time.unit = units$time, value.unit = units$value,
                   xlab = xlab, pretty.x.breaks = pretty_x_breaks,
                   poly.alpha = poly.alpha, cex = cex, sim.lwd = sim.lwd, error.lwd =error.lwd,
                   add.legend.text = add_led,
                   legend.plot.args = legend_plot_args,
                   main.plot.args = main_plot_args, 
                   ylab = NA, show.main = T, show.legend = T, 
                   rm.zero.neg.rows = F, ymax.rel.add = 0.0,
                   legend_format = legend_format,
                   panel.first = panel_fn(profile, counter, is_fraction, "linear"), ...)
      
      dev.off()
    }
    
    # log plots
    p_log <- plot_log_fn(profile, counter, is_fraction)
    msg_fn(paste("\tLog:", p_log), nl = FALSE)
    if (p_log) {
      file_name <- paste0(id_fn(profile, counter, is_fraction, "log"), ".", format)
      full_name <- file.path(plot_folder, file_name)
      
      plot_fn(full_name)
      par_plot_fn()
      
      plot.matched(profile, observed_data,
                   avg.fn = avg.fn, min.var.fn = min.var.fn, max.var.fn = max.var.fn,
                   time.unit = units$time, value.unit = units$value,
                   xlab = xlab, pretty.x.breaks = pretty_x_breaks,
                   poly.alpha = poly.alpha, cex = cex, sim.lwd = sim.lwd, error.lwd =error.lwd,
                   add.legend.text = add_led,
                   legend.plot.args = legend_plot_args,
                   main.plot.args = append(main_plot_args,  list(log = "y")), 
                   ylab = NA, show.main = T, show.legend = T, 
                   rm.zero.neg.rows = F, ymax.rel.add = 0.0,
                   legend_format = legend_format,
                   panel.first = panel_fn(profile, counter, is_fraction, "log"), ...)
      
      dev.off()
    }
    msg_fn("")
  } 
  
  msg_fn(paste("Plotted <", length(profiles), "> profiles"))
}
