
library(sitar)
library(tidyverse)

get_pBMI <- function(df, age_limit) {
  tmp <- df %>% filter(age <= age_limit) %>% split(.$id)
  map_df(tmp, ~ getPeak(x = .x$age, y = .x$cbmi))
}

get_rBMI <- function(df, age_min) {
  tmp <- df %>% filter(age > age_min) %>% split(.$id)
  map_df(tmp, ~ getTrough(x = .x$age, y = .x$cbmi))
}

get_pv <- function(
    data,
    id_col = "id",
    age_col = "age",
    value_col = "height",
    n_window = 1,
    spar = NULL) {
  
  individual_list <- split(data, data[[id_col]])
  
  results_list <- lapply(individual_list, function(dat) {
    x <- dat[[age_col]]
    y <- dat[[value_col]]
    
    if (length(unique(x)) < 4L) {
      return(data.frame(
        id = unique(dat[[id_col]]),
        APV = NA_real_,
        PV = NA_real_
      ))
    }
    
    ss <- smooth.spline(x, y, spar = spar)
    
    vfun <- function(a) as.numeric(predict(ss, a, deriv = 1L)$y)
    
    # Find global maximum velocity (allowing boundary)
    rng <- range(x, na.rm = TRUE)
    opt <- optimize(vfun, lower = rng[1], upper = rng[2], maximum = TRUE)
    peak_age_guess <- opt$maximum
    peak_vel_guess <- opt$objective
    
    # Quadratic refinement
    vel_all <- vfun(sort(unique(x)))
    xy <- data.frame(x = sort(unique(x)), y = vel_all)
    idx <- which.min(abs(xy$x - peak_age_guess))
    low <- max(1, idx - n_window)
    high <- min(nrow(xy), idx + n_window)
    fit <- lm(y ~ poly(x, 2, raw = TRUE), data = xy[low:high, ])
    
    if (fit$rank < 3) {
      peak_age_refined <- peak_age_guess
      peak_vel_refined <- peak_vel_guess
    } else {
      peak_age_refined <- -coef(fit)[2] / (2 * coef(fit)[3])
      peak_vel_refined <- predict(fit, newdata = data.frame(x = peak_age_refined))
      # Safety: if refined peak is outside observed range, revert to guess
      if (peak_age_refined < min(x) || peak_age_refined > max(x)) {
        peak_age_refined <- peak_age_guess
        peak_vel_refined <- peak_vel_guess
      }
    }
    
    data.frame(
      id = unique(dat[[id_col]]),
      APV = peak_age_refined,
      PV = peak_vel_refined
    )
  })
  
  do.call(rbind, results_list)
}