# helper
.geo.sd <- function(x, na.rm = TRUE) {
  x[x <= 0] <- 1E-99
  exp(stats::sd(log(x), na.rm = na.rm))
}

# variance functions
q16 <- function(values) stats::quantile(values, probs = 0.16, names = F)
q84 <- function(values) stats::quantile(values, probs = 0.84, names = F)

std.dev.min <- function(values) mean(values) - stats::sd(values)
std.dev.max <- function(values) mean(values) + stats::sd(values)

geo.dev.min <- function(values) geo.mean(values) * .geo.sd(values)
geo.dev.max <- function(values) geo.mean(values) / .geo.sd(values)

# average functions
geo.mean <- function(x, na.rm = TRUE) {
  x[x <= 0] <- 1E-99
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
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
range.profile <- function(profile, range.type = c("time", "data")) {
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

  min <- NA
  max <- NA
  min <- min(data, na.rm = TRUE)
  max <- max(data, na.rm = TRUE)

  return(c(min, max))
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
trim.time <- function(profile, from = NA, to = NA) {
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
  idx <- which(times >= from & times <= to)
  idx <- sort(idx)

  result <- profile
  result$data <- result$data[idx, ]

  return(result)
}

merge.profile <- function(from, into, meta.data = c("into, from"),
                          units = c("check, into, from"),
                          keep.duplicates = c("into", "from"),
                          validate.result = T) {

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
    from.duplicates <- which(from$data$Time %in% into$data$Time)
    if (length(from.duplicates) > 0)
      from$data <- from$data[-from.duplicates,]
  } else {
    into.duplicates <- which(into$data$Time %in% from$data$Time)
    if (length(into.duplicates) > 0)
      into$data <- into$data[-into.duplicates,]
  }

  # bind data
  res.data <- rbind(into$data, from$data)
  res.data <- res.data[order(res.data$Time),]

  result$data <- res.data
  # check for valid results
  if (validate.result) {
    if (!is.valid(result))
      stop("Result is not a valid profile")
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
  patter.range <- range.profile(in.profile, range.type = "time")
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
    result.profile <- merge(result.profile, in.profile)
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


calculate.all.pred.obs <- function(matched.list, obs.data, time.unit, value.unit,
                                   only.obs.times = F,
                                   interpol.method = c("linear", "spline"),
                                   interpol.spline.method = c("fmm", "periodic", "natural", "monoH.FC", "hyman"),
                                   auc.method = c("linear", "linlog", "spline")) {

  auc.method <- match.arg(auc.method)

  result <- data.frame()
  for (match in matched.list) {
    message(paste("Processing id <", match$id, ">"))
    tmp <- calculate.pred.obs(match, obs.data, time.unit, value.unit,
                              only.obs.times = only.obs.times,
                              interpol.method = interpol.method,
                              interpol.spline.method = interpol.spline.method,
                              auc.method = auc.method)
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
                               auc.method = c("linear", "linlog", "spline")) {

  auc.method <- match.arg(auc.method)

  obs <- .gather.ids(obs.data, matched$obs.ids)
  if (length(obs) == 0) {
    stop(paste("Attention: Matched profiles with id < ", matched$id, "> do not have any observed data"))
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

    obs.max <- calculate.max(od)
    pred.max <- calculate.max(pro)

    obs.range <- range(od, range.type = "time")
    obs.auc <- calculate.auc(od, type = auc.method)
    pred.auc <- calculate.auc(trim.time(pro, from = obs.range[1], to = obs.range[2]),
                              type = auc.method)

    tmp.res <- data.frame(
      sim.id = matched$id,
      exp.id = od$id,
      group = od$group,
      ref = od$reference,
      mol.name = pro.mol$name,
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

cor.sspd <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  mdlq <- cor.mdlq(pred, obs)
  sign(mdlq) * (exp(mdlq) - 1)
}

cor.zeta <- function(pred, obs) {
  pred <- .rm.units(pred)
  obs <- .rm.units(obs)
  exp(stats::median(abs(log(pred / obs))) - 1)
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
      "SSPB" = units::set_units(100 * cor.sspd(pred, obs), "%"), # in %
      "Zeta" = units::set_units(100 * cor.zeta(pred, obs), "%") # in %
    )
    results[[metric.name]] <- tmp
  }

  return(results)
}
