
# get n distinct colors
.get_distinct_colors <- function(n) {
  set.seed(8675309)
  return(randomcoloR::distinctColorPalette(n))
}

# get n distinct pchs
.get_distinct_pch <- function(n, show.warnings = TRUE) {
  
  # acceptable pchs
  # solids (4) + hallow (4) + crosses (4)
  acceptable <- c(19, 17, 15, 18, 1, 2, 0, 5, 13, 6, 12, 9)
  
  if (n > length(acceptable))
    warning("Number of distinct pchs is too large. Point types will be recycled.")
  
  return(rep(acceptable, length.out = n))
}


plot_gof_pk <- function(pred.obs.data, #
                        pk.parameter, #
                        lab.text, #
                        group_by = "ref", #
                        subgroup_by = NA, #
                        unit.lab = NULL, # NA for no unit
                        col = NA, #
                        pch = NA, #
                        min = NA, #
                        max = NA, #
                        nice.min = T, #
                        nice.max = T, #
                        symmetry = F, #
                        lwd = 2.5, #
                        legend.cex = 1.25, #
                        legend.ncol = 1, #
                        spit.legend.at = NULL, #
                        legend.titles = NULL, #
                        show.loi = TRUE, #
                        show.2.fold = TRUE, #
                        show.1.25.fold = TRUE, #
                        show.guest = FALSE, #
                        guest.delta = 1, #
                        par_fn = NULL,
                        format = c("dev", "pdf", "tiff", "png", "jpeg"),
                        format_opts = NULL,
                        out_file = NULL,
                        ...) {
  
  # graphics reset
  on.exit(graphics::layout(1))
  opar<-par(no.readonly = TRUE)
  on.exit(suppressWarnings(par(opar)), add = TRUE, after=FALSE)
  
  # inits
  par_plot_fn <- function() {
    
    if (!is.null(par_fn))
      par_fn()
    else
      par(mar = c(5.1, 4.1, 4.1, 2.1) + c(0, 1, -1.5, 0))
  }
  
  format <- match.arg(format)
  plot_fn <- function() {
    if (format != "dev")
      do.call(format, append(list(out_file), format_opts))
  }
  
  plot_fn_end <- function() {
    if (format != "dev")
      dev.off()
  }
  
  
  ################
  pk.parameter <- tolower(pk.parameter)
  target.pred <- paste0("pred.", pk.parameter)
  target.obs <- paste0("obs.", pk.parameter)

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
                                            b = lab.text,
                                            c = lab_unit), "[[", 1))
  
  xlab <- substitute(a ~ b ~ c, lapply(list(a = "Observed",
                                            b = lab.text,
                                            c = lab_unit), "[[", 1))
  plot_fn()
  par_plot_fn()
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
  if (show.loi) {
    graphics::abline(0, 1, lwd = lwd)
  }
    
  if (show.2.fold) {
    graphics::abline(0 , 2, lty = 2, untf = T, lwd = lwd)
    graphics::abline(0 , 0.5, lty = 2, untf = T, lwd = lwd)
  }
  
  if (show.1.25.fold) {
    graphics::abline(0, 1.25, lty = 3, untf = T, lwd = lwd)
    graphics::abline(0, 1/1.25, lty = 3, untf = T, lwd = lwd)
  }
  
  if (show.guest) {
    graphics::curve(guest.delta + 2*(x -1), add = TRUE, from= 1, to = range[2] * 10, lwd=lwd)
    graphics::curve(1/(guest.delta + 2*(1/x -1)), add = TRUE, to= 1, from = range[1] / 10, lwd=lwd)
    graphics::curve(x**2 / (guest.delta + 2*(x -1)), add = TRUE, from= 1, to = range[2] * 10, lwd=lwd)
    graphics::curve(1/((1/x)**2 / (guest.delta + 2*(1/x -1))), add = TRUE, to = 1, 
                    from = range[1] / 10, lwd=lwd)
  }
  
  df <- cbind(df, pred.obs.data[group_by])
  colnames(df) <- c("pred", "obs", "group")
  if (!is.na(subgroup_by)) {
    df <- cbind(df, pred.obs.data[subgroup_by])
    colnames(df) <- c("pred", "obs", "group", "subgroup")
  }
  
  groups <- unique(df$group)
  n.groups <- length(groups)
  
  col_df <- col
  if (is.data.frame(col_df) || (length(col) <= 1 && is.na(col)))
    col <- .get_distinct_colors(n.groups)
  
  if (length(col) < n.groups)
    col <- rep(col, length.out = n.groups)
  
  
  pch.n <- n.groups
  if (!is.na(subgroup_by))
    pch.n <- max(pch.n, length(unique(df[["subgroup"]])))
  
  if (length(pch) <= 1 && is.na(pch))
    pch <- .get_distinct_pch(pch.n)
  
  if (length(pch) < pch.n)
    col <- rep(pch, length.out = pch.n)
  
  cex <- unlist(list(...)["cex"])
  cex <- if (is.null(cex)) 1.0 else cex
  
  
  if (is.data.frame(col_df)) {
    col <- .extract_color_map_col(col_df, groups)
  }
  
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
          new_color <- new_col[i] 
          
          
          pch <- c(pch, pch_map[[sub]])
          col <- c(col, new_color)
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
  
  if (!is.null(spit.legend.at)) {
    if (is.numeric(spit.legend.at)) {
      sp_idx <- spit.legend.at:length(groups)
    } else {
      sp_idx <- grep(spit.legend.at, groups)
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
  
  plot_fn_end()
}


color_map <- function(df, ref_column = "ref", col = NULL) {
  
  if (is.list(df))
    df <- df$data
  
  reference <- df[[ref_column]]
  if (is.null(reference))
    stop("ref_column is not a column of df", call. = FALSE)
    
  reference <- unique(reference)
  n <- length(reference)
  
  if (is.null(col))
    col <- .get_distinct_colors(n)
  
  if (length(col) < n)
    stop(paste("col must be of length <", n, ">"), call. = FALSE)
  
  
  result <- data.frame(reference = reference, col = col[1:n])
  return(result)
}

.extract_color_map_col <- function(map, keys) {
  cols <- c()
  for (key in keys) {
    cols <- c(cols, map %>% dplyr::filter(reference == key) %>% pull(col))
  }
  return(cols)
}

