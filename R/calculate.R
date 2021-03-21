# helper
.geo.sd <- function(x, na.rm = TRUE) {

  if (any(x <= 0)) {
    x <- x + 1
  }
  mean <- geo.mean(x, na.rm = na.rm)

  exp(sqrt(sum(log(x/mean)**2) / length(x)))
}

# variance functions
q16 <- function(values) stats::quantile(values, probs = 0.16, names = F)
q84 <- function(values) stats::quantile(values, probs = 0.84, names = F)

q05 <- function(values) stats::quantile(values, probs = 0.05, names = F)
q95 <- function(values) stats::quantile(values, probs = 0.95, names = F)

std.dev.min <- function(values) mean(values) - stats::sd(values)
std.dev.max <- function(values) mean(values) + stats::sd(values)

geo.dev.min <- function(values) geo.mean(values) * .geo.sd(values)
geo.dev.max <- function(values) geo.mean(values) / .geo.sd(values)

# average functions
geo.mean <- function(x, na.rm = TRUE) {
  has_zeros <- FALSE
  if (any(x <= 0)) {
    x <- x + 1
    has_zeros = TRUE
  }

  result <- exp(sum(log(x), na.rm = na.rm) / length(x))
  if (has_zeros)
    result <- result - 1

  return(result)
}

# individual population profiles -> mean population profiles
average.pop.profile <- function(profile,
                                avg.fn = base::mean,
                                min.var.fn = std.dev.min,
                                max.var.fn = std.dev.max) {
  if (!is.profile(profile)) {
    stop("Input must be of class profile")
  }

  if (profile$data.type != "individual") {
    stop("Input must be of data.type individual")
  }

  if (profile$type != "population") {
    stop("Input must be of type population")
  }

  result <- profile

  result$data <- data.frame(Time = profile$data$Time)
  data <- profile$data[-1]

  avg <- apply(data, 1, avg.fn)
  min <- apply(data, 1, min.var.fn)
  max <- apply(data, 1, max.var.fn)

  result$data <- data.frame(
    Time = profile$data$Time,
    Avg = avg,
    Min = min,
    Max = max
  )

  result$data.type <- "mean"
  return(result)
}

# calculate ranges of a profile
range <- function(...) UseMethod("range")
range.profile <- function(profile, range.type = c("time", "data"),
                          smart.md = FALSE,
                          smart.md.threshold = 15) {
  if (!is.profile(profile)) {
    stop("Input must be of class profile")
  }

  type <- match.arg(range.type)

  data <- NA
  if (type == "time") {
    data <- profile$data$Time
  } else {
    data <- profile$data[-1]
  }

  if (type == "time" && smart.md && length(data > 3)) {
    d <- diff(data)
    m_d <- mean(d)
    groups <- which(d >= m_d * smart.md.threshold)
    if (length(groups) > 0) {
      groups <- c(0, groups, length(data))
      if (any(diff(groups) == 1))
      {
        message(paste("Smart MD found <", length(groups) - 1 , "> groups",
        "but with at least one single point group. Range will be computed without grouping."))
      }
      else {
      message(paste("Smart MD found <", length(groups) - 1 , "> groups"))

      df <- data.frame(from = groups[-length(groups)] + 1, to = groups[-1])
      res <- data.frame(t(apply(df, 1, function(x) c(data[x[1]], data[x[2]]))))
      colnames(res) <- c("from", "to")
      return(res)
      }
    }
  }

  return(base::range(data, na.rm = TRUE))
}

convert <- function(...) UseMethod("convert")
convert.profile <- function(profile, time.unit = NA, value.unit = NA) {
  if (!is.profile(profile)) {
    stop("Input must be of class profile")
  }

  result <- profile
  if (!is.na(time.unit)) {
    if (is.character(time.unit)) {
      time.unit <- units::as_units(time.unit)
    }

    if (is.null(result$time.unit)) {
      stop("No time unit set in profile")
    }

    result$data$Time <- .convert.units(
      result$data$Time,
      result$time.unit, time.unit
    )
    result$time.unit <- time.unit
  }

  if (!is.na(value.unit)) {
    if (is.character(value.unit)) {
      value.unit <- units::as_units(value.unit)
    }

    if (is.null(result$value.unit)) {
      stop("No value unit set in profile")
    }

    result$data[-1] <- .convert.units(
      result$data[-1],
      result$value.unit, value.unit, profile$molecule$MW
    )
    result$value.unit <- value.unit
  }

  return(result)
}

# trims a data.set for a certain time interval
trim.time <- function(profile, from = NA, to = NA, tol= .Machine$double.eps^0.5) {
  if (!is.profile(profile)) {
    stop("Input must be of class profile")
  }


  if (is.na(from)) {
    from <- -Inf
  }

  if (is.na(to)) {
    to <- +Inf
  }

  times <- profile$data$Time
  idx <- which(times >= from - tol & times <= to + tol)
  idx <- sort(idx)

  result <- profile
  result$data <- result$data[idx, ]

  return(result)
}

.test.near <- function(a,b) {
  dplyr::near(a,b, .Machine$double.eps**(0.25))
}

.find.near.duplicates <- function(a, in_b) {
  which(unlist(sapply(a, function(x) {
    matches <- .test.near(x, in_b)
    any(matches)
  })))
}

merge.profile <- function(from, into, meta.data = c("into, from"),
                          units = c("check, into, from"),
                          keep.duplicates = c("into", "from"),
                          validate.result = T) {

  from$data$Time <- from$data$Time
  into$data$Time <- into$data$Time

  if (!is.profile(from))
    stop("from must be of class profile")

  if (!is.profile(into))
    stop("into must be of class profile")

  if (ncol(from$data) != ncol(into$data))
    stop("from and into are not compatible (data column count)")

  meta.data <- match.arg(meta.data)
  units <- match.arg(units)
  keep.duplicates <- match.arg(keep.duplicates)

  # convert units
  if (units == "check") {
    if (from$time.unit != into$time.unit || from$value.unit != into$value.unit) {
      stop("Units are not equal")
    }
  } else if (units == "into") {
    from <- convert.profile(from, time.unit = into$time.unit, value.unit = into$value.unit)
  } else if (units == "from") {
    into <- convert.profile(into, time.unit = from$time.unit, value.unit = from$value.unit)
  }

  # apply meta.data
  result <- if (meta.data == "into") into else from

  # handle duplicates
  if (keep.duplicates == "into") {
    #from.duplicates <- which(from$data$Time %in% into$data$Time)
    from.duplicates <- .find.near.duplicates(from$data$Time, into$data$Time)
    if (length(from.duplicates) > 0)
      from$data <- from$data[-from.duplicates,]
  } else {
    #into.duplicates <- which(into$data$Time %in% from$data$Time)
    into.duplicates <- .find.near.duplicates(into$data$Time, from$data$Time)
    if (length(into.duplicates) > 0)
      into$data <- into$data[-into.duplicates,]
  }

  # bind data
  res.data <- rbind(into$data, from$data)
  res.data <- res.data[order(res.data$Time),]
  result$data <- res.data

  # check for valid results
  if (validate.result) {
    if (!is.valid(result, msg = "message")) {
      stop("Result is not a valid profile")
    }
  }

  return(result)
}


interpol.profile <- function(in.profile,
                             pattern.profile,
                             method = c("linear", "spline"),
                             spline.method = c("fmm", "periodic", "natural", "monoH.FC", "hyman"),
                             conserve.time.unit = F,
                             only.pattern.times = F) {
  method <- match.arg(method)
  spline.method <- match.arg(spline.method)

  if (!is.profile(in.profile))
    stop("in.profile must be of class profile")

  if (!is.profile(pattern.profile))
    stop("pattern.profile must be of class profile")

  # convert time
  old.time.unit <- in.profile$time.unit
  in.profile <- convert.profile(in.profile, time.unit = pattern.profile$time.unit)

  # check for time range
  in.time.range <- range.profile(in.profile, range.type = "time")
  patter.range <- range.profile(pattern.profile, range.type = "time")
  if (patter.range[1] < in.time.range[1] && .test.near(in.time.range[1], patter.range[1])) {
    message("Profile Min times larger than Obs Min but close. Will try to fix this.")
    in.profile$data$Time[0] = patter.range[1]
  }
  if (patter.range[2] > in.time.range[2] && .test.near(in.time.range[2], patter.range[2])) {
    message("Profile Max times less than Obs Max but close. Will try to fix this.")
    in.profile$data$Time[length(in.profile$data$Time)] = patter.range[2]
  }

  in.time.range <- range.profile(in.profile, range.type = "time")
  patter.range <- range.profile(pattern.profile, range.type = "time")
  if (patter.range[1] < in.time.range[1] || patter.range[2] > in.time.range[2]) {
    stop("Pattern time range is not inside of in.profile range")
  }

  # get times
  pattern.times <- pattern.profile$data$Time
  in.times <- in.profile$data$Time

  result.profile <- in.profile
  res.data <- data.frame(Time = pattern.times)
  for (col in 2:ncol(in.profile$data)) {
    in.col <- in.profile$data[, col]
    if (is.na(in.col[1])) {
      res.data <- cbind(res.data, NA)
    } else {
      if (method == "linear") {
        approx.data <- stats::approx(in.times, in.col, pattern.times)$y
      } else {
        spline.fn <- stats::splinefun(in.times, in.col, method = spline.method)
        approx.data <- spline.fn(pattern.times)
      }
      res.data <- cbind(res.data, approx.data)
    }
  }

  colnames(res.data) <- colnames(result.profile$data)
  result.profile$data <- res.data

  if (!only.pattern.times) {
    result.profile <- merge.profile(result.profile, in.profile)
  }

  if (conserve.time.unit) {
    result.profile <- convert.profile(result.profile, time.unit = old.time.unit)
  }

  return(result.profile)
}


.calculate_auc_linlog <- function(times, values) {

  if (length(times) != length(values))
    stop("Error: times and values must have the same length")

  if (length(times) < 2)
    return(0.)

  auc = 0.

  for (i in 2:length(times)) {
    t1 <- times[i - 1]
    t2 <- times[i]
    if (dplyr::near(t1, t2)) {
      next
    }

    c1 <- values[i - 1]
    c2 <- values[i]
    if (c1 > c2 && c2 > 0.)
      auc <- auc + (c1 - c2) * (t2 - t1)/(log(c1) - log(c2))
    else
      auc <- auc + 0.5 * (c1 + c2) * (t2 - t1)
  }



  return(auc)
}

calculate.auc <- function(profile, type = c("linear", "linlog", "spline")) {
  type <- match.arg(type)
  if (!is.profile(profile)) {
    stop("Input must be of class profile")
  }

  if (profile$data.type != "mean") {
    stop("Input must be of data.type mean")
  }

  time.unit <- profile$time.unit
  value.unit <- profile$value.unit

  if (type == "linear" || type == "spline") {
    auc <- MESS::auc(profile$data$Time, profile$data$Avg, type = type)
  }
  else {
    auc <- .calculate_auc_linlog(profile$data$Time, profile$data$Avg)
  }

  units(auc) <- value.unit * time.unit
  return(auc)
}

calculate.max <- function(profile) {
  if (!is.profile(profile)) {
    stop("Input must be of class profile")
  }

  if (profile$data.type != "mean") {
    stop("Input must be of data.type mean")
  }

  time.unit <- profile$time.unit
  value.unit <- profile$value.unit

  max.col <- which(profile$data$Avg == max(profile$data$Avg, na.rm = TRUE))[1]

  t.max <- profile$data$Time[max.col]
  v.max <- profile$data$Avg[max.col]

  units(t.max) <- time.unit
  units(v.max) <- value.unit


  return(list(t.max = t.max, value.max = v.max))
}

all_pred_vs_obs <- function(matched.list, obs.data,
                            time.unit, value.unit,
                            interpol.method = c("linear", "spline"),
                            interpol.spline.method = c("fmm", "periodic", "natural", "monoH.FC", "hyman"),
                            auc.method = c("linear", "linlog", "spline")) {

  results <- list(meta = list(time.unit = time.unit, value.unit = value.unit),
                  data = data.frame())

  for (match in matched.list) {
    message(paste("Processing id <", match$id, ">"))
    tmp <- pred_vs_obs(match, obs.data, time.unit, value.unit,
                              interpol.method = interpol.method,
                              interpol.spline.method = interpol.spline.method)
    results$data <- rbind(results$data, tmp$data)
  }

  return(results)
}


pred_vs_obs <- function(matched, obs.data,
                        time.unit,
                        value.unit,
                        interpol.method = c("linear", "spline"),
                        interpol.spline.method = c("fmm", "periodic", "natural", "monoH.FC", "hyman")) {

  obs <- .gather.ids(obs.data, matched$obs.ids)

  results <- list(meta = list(time.unit = time.unit, value.unit = value.unit),
                  data = data.frame())
  for (pro in matched$profiles) {
    pro.mol <- pro$molecule
    match <- NA
    for (od in obs) {
      if (od$molecule$id == pro.mol$id) {
        match <- od
        break
      }
    }

    if (!is.profile(match)) {
      message(paste("Unexpected molecule mismatch in < ", matched$id, ">"))
      next
    }

    pro <- convert.profile(pro, time.unit = time.unit, value.unit = value.unit)
    od <- convert.profile(od, time.unit = time.unit, value.unit = value.unit)

    pro <- interpol.profile(pro, od,
                            method = interpol.method,
                            spline.method = interpol.spline.method,
                            only.pattern.times = T)

    result_df <- data.frame(Time = pro$data$Time,
                            Pred = pro$data$Avg,
                            Obs = od$data$Avg,
                            sim.id = matched$id,
                            exp.id = od$id,
                            ref = od$reference,
                            cite.key = od$citekey,
                            dose = od$dose,
                            dose.unit = od$dose.unit,
                            admin.route = od$route,
                            group = od$group,
                            group2 = od$group2,
                            group3 = od$group3,
                            ref = od$reference,
                            mol.name = pro.mol$name,
                            mol.display.name = pro.mol$display.name,
                            mol.id = pro.mol$id)

    results$data <- rbind(results$data, result_df)
  }
  return(results)
}



calculate.all.pred.obs <- function(matched.list, obs.data, time.unit, value.unit,
                                   only.obs.times = F,
                                   interpol.method = c("linear", "spline"),
                                   interpol.spline.method = c("fmm", "periodic", "natural", "monoH.FC", "hyman"),
                                   auc.method = c("linear", "linlog", "spline"),
                                   smart.md.threshold = 15) {

  auc.method <- match.arg(auc.method)

  result <- data.frame()
  for (match in matched.list) {
    message(paste("Processing id <", match$id, ">"))
    tmp <- calculate.pred.obs(match, obs.data, time.unit, value.unit,
                              only.obs.times = only.obs.times,
                              interpol.method = interpol.method,
                              interpol.spline.method = interpol.spline.method,
                              auc.method = auc.method,
                              smart.md.threshold = smart.md.threshold)
    result <- rbind(result, tmp)
  }

  return(result)
}


calculate.pred.obs <- function(matched, obs.data,
                               time.unit,
                               value.unit,
                               only.obs.times = F,
                               interpol.method = c("linear", "spline"),
                               interpol.spline.method = c("fmm", "periodic", "natural", "monoH.FC", "hyman"),
                               auc.method = c("linear", "linlog", "spline"),
                               smart.md.threshold = 15) {

  auc.method <- match.arg(auc.method)

  obs <- .gather.ids(obs.data, matched$obs.ids)
  if (length(obs) == 0) {
    stop(paste("Attention: Matched profiles with id < ", matched$id, "> do not have any observed data"))
  }

  if (is.null(obs[[1]])) {
    stop(paste("Attention: Matched profiles with id < ", matched$id, "> could not match observed data."))
  }

  if (!is.na(matched$obs.ids) && length(obs) != length(matched$obs.ids)) {
    stop(paste("Attention: Matched profiles with id < ", matched$id, "> have missing observed data"))
  }

  results <- data.frame()
  for (pro in matched$profiles) {
    pro.mol <- pro$molecule
    match <- NA
    for (od in obs) {
      if (od$molecule$id == pro.mol$id) {
        match <- od
        break
      }
    }

    if (!is.profile(match)) {
      message(paste("Unexpected molecule mismatch in < ", matched$id, ">"))
      next
    }
    pro <- convert.profile(pro, time.unit = time.unit, value.unit = value.unit)
    od <- convert.profile(od, time.unit = time.unit, value.unit = value.unit)

    pro <- interpol.profile(pro, od,
                            method = interpol.method,
                            spline.method = interpol.spline.method,
                            only.pattern.times = only.obs.times)


    # max values
    obs.range_full <- range(od, range.type = "time", smart.md = FALSE,
                       smart.md.threshold = smart.md.threshold)
    pro_trimmed <- trim.time(pro, from = obs.range_full[1], to = obs.range_full[2])
    obs.max <- calculate.max(od)
    pred.max <- calculate.max(pro_trimmed)

    # AUC
    obs.range <- range(od, range.type = "time", smart.md = TRUE,
                       smart.md.threshold = smart.md.threshold)
    if (!is.null(nrow(obs.range))) {
      split.obs <- apply(obs.range, 1, function(x) trim.time(od, x[1], x[2]))
      split.pro <- apply(obs.range, 1, function(x) trim.time(pro, x[1], x[2]))

      obs.auc <- sum(sapply(split.obs, calculate.auc, type = auc.method))
      pred.auc <- sum(sapply(split.pro, calculate.auc, type = auc.method))

    } else {
      obs.auc <- calculate.auc(od, type = auc.method)
      pred.auc <- calculate.auc(pro_trimmed, type = auc.method)
    }

    tmp.res <- data.frame(
      sim.id = matched$id,
      exp.id = od$id,
      group = od$group,
      group2 = od$group2,
      group3 = od$group3,
      ref = od$reference,
      cite.key = od$citekey,
      dose = od$dose,
      dose.unit = od$dose.unit,
      admin.route = od$route,
      mol.name = pro.mol$name,
      mol.display.name = pro.mol$display.name,
      mol.id = pro.mol$id,
      obs.value.max = obs.max$value.max,
      pred.value.max = pred.max$value.max,
      obs.t.max = obs.max$t.max,
      pred.t.max = pred.max$t.max,
      obs.auc = obs.auc,
      pred.auc = pred.auc
    )

    results <- rbind(results, tmp.res)
  }

  return(results)
}


cor.mse <- function(pred, obs) sum((pred - obs)**2) / length(pred)

cor.rmse <- function(pred, obs) sqrt(cor.mse(pred, obs))

cor.mae <- function(pred, obs) sum(abs(pred - obs)) / length(pred)

cor.mdae <- function(pred, obs) stats::median(abs(pred - obs))

cor.mape <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  sum(abs((pred - obs) / obs)) / length(pred)
}

cor.smape <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  sum(abs((pred - obs) / (0.5 * (pred + obs)))) / length(pred)
}

cor.me <- function(pred, obs) sum(pred - obs) / length(pred)

cor.mpe <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  sum((pred - obs) / obs) / length(pred)
}

cor.mdlq <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  stats::median(log(pred / obs))
}

cor.sspb <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  mdlq <- cor.mdlq(pred, obs)
  sign(mdlq) * (exp(abs(mdlq)) - 1)
}

cor.zeta <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  exp(stats::median(abs(log(pred / obs)))) - 1
}

cor.max_ae <- function(pred, obs) max(abs(pred - obs))
cor.max_ape <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  max(abs(pred - obs) / obs)
}

cor.min_ae <- function(pred, obs) min(abs(pred - obs))

cor.min_ape <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  min(abs(pred - obs) / obs)
}


cor.mrd <- function(pred, obs) {
  pred_u <- .rm.units(pred)
  obs_u <- .rm.units(obs)

  res <- 10**(sqrt(sum((log10(pred_u) - log10(obs_u))**2) / length(pred_u)))
  if (.has.units(pred) && .has.units(obs)) {
    units(res) <- units(pred)
  }
  return(res)
}

cor.gmfe <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)

  10**(sum(abs(log10(pred / obs))) / length(pred))
}

cor.metric.table <- function(metrics) {
  stats::setNames(data.frame(matrix(unlist(metrics),
    nrow = length(metrics), byrow = T
  ),
  row.names = names(metrics)
  ), names(metrics[[1]]))
}

cor.metrics.pred.obs <- function(df, groups = NULL) {
  if (length(groups) > 0) {
    df <- df[which(df$group %in% groups), ]
  }

  if (nrow(df) < 2) {
    stop("Input must have at least 2 entry rows")
  }

  names <- colnames(df)
  obs.idx <- which(grepl("obs.", names))

  results <- list()
  for (idx in obs.idx) {
    obs <- df[, idx]
    pred <- df[, idx + 1]

    metric.name <- sub("obs.", "", names[idx])

    tmp <- data.frame(
      "MRD" = cor.mrd(pred, obs),
      "GMFE" = cor.gmfe(pred, obs),
      "Max AE" = cor.max_ae(pred, obs),
      "Min AE" = cor.min_ae(pred, obs),
      "Max APE" = units::set_units(100 * cor.max_ape(pred, obs), "%"), # in %
      "Min APE" = units::set_units(100 * cor.min_ape(pred, obs), "%"), # in %
      "MSE" = cor.mse(pred, obs),
      "RMSE" = cor.rmse(pred, obs),
      "MAE" = cor.mae(pred, obs),
      "MdAE" = cor.mdae(pred, obs),
      "MAPE" = units::set_units(100 * cor.mape(pred, obs), "%"), # in %
      "sMAPE" = units::set_units(100 * cor.smape(pred, obs), "%"), # in %
      "ME" = cor.me(pred, obs),
      "MPE" = units::set_units(100 * cor.mpe(pred, obs), "%"), # in %
      "MdLQ" = cor.mdlq(pred, obs),
      "SSPB" = units::set_units(100 * cor.sspb(pred, obs), "%"), # in %
      "Zeta" = units::set_units(100 * cor.zeta(pred, obs), "%") # in %
    )
    results[[metric.name]] <- tmp
  }

  return(results)
}



calculate_ratios <- function(df, id.group = "group", effect.group = "group2", 
                             control.name = "control", effect.name = "ddi",
                             ref = c("control", "effect")) {
  
  ref <- match.arg(ref)
  
  ids <- unique(df[[id.group]])
  effects <- unique(df[[effect.group]])
  
  # basic error handling
  if (is.null(ids))
    stop(paste("Could not find id.group <", id.group, "> in data.frame"), call. = FALSE)
  
  if (is.null(effects))
    stop(paste("Could not find effect.group <", effect.group, "> in data.frame"), call. = FALSE)
  
  if (length(effects) != 2)
    stop("Number of unique effect names in effect.group must be 2", call. = FALSE)
  
  if (!(control.name %in% effects))
    stop(paste("Did not find control.name <", control.name, "> in effect.group"), call. = FALSE)
  
  if (!(effect.name %in% effects))
    stop(paste("Did not find effect.name <", effect.name, "> in effect.group"), call. = FALSE)
  
  result <- data.frame()
  for (id in ids) {
    tmp <- df[df[[id.group]] == id, ]
    
    effect.row <- tmp[tmp[[effect.group]] == effect.name, ]
    control.row <- tmp[tmp[[effect.group]] == control.name, ]
    
    if(nrow(control.row) != 1)
      stop(paste("Error for id.group <", id, ">: Expected one control entry"), call. = FALSE)
    
    if(nrow(effect.row) < 1)
      stop(paste("Error for id.group <", id, ">: Expected at least one effect entry"), call. = FALSE)
    
    n_eff_rows <- nrow(effect.row)
    for (i in 1:n_eff_rows) {
    
      effect <- effect.row[i,]
      row <- if(ref == "control") control.row else effect
      
      row$obs.value.max  <- effect$obs.value.max / control.row$obs.value.max
      row$pred.value.max  <- effect$pred.value.max / control.row$pred.value.max
      row$obs.t.max  <- effect$obs.t.max / control.row$obs.t.max
      row$pred.t.max  <- effect$pred.t.max / control.row$pred.t.max
      row$obs.auc  <- effect$obs.auc / control.row$obs.auc
      row$pred.auc  <- effect$pred.auc / control.row$pred.auc
      
      result <- rbind(result, row)
    }
  }
  
  return(result)
}
