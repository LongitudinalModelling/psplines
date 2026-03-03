
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

get_phv <- function(df, age_limit) {
  tmp <- df %>% filter(age <= age_limit) %>% split(.$id)
  map_df(tmp, ~ getPeak(x = .x$age, y = .x$cht, Dy = T))
}

get_pwv <- function(df, age_limit) {
  tmp <- df %>% filter(age <= age_limit) %>% split(.$id)
  map_df(tmp, ~ getPeak(x = .x$age, y = .x$cwt, Dy = T))
}

get_vel <- function(df, var) {
  df %>% select(id, age, !!sym(var)) %>% nest_by(id) %>% mutate(
    deriv = list({
      x <- data$age
      y <- data[[var]]
      ss <- smooth.spline(x, y)
      ages_m <- c(1, 6, 12, 24)
      pred <- predict(ss, x = sqrt(ages_m / 12), deriv = 1L)
      tibble(age_m = ages_m, !!paste0("v", var) := pred$y)})) %>%
    select(id, deriv) %>% unnest(cols = deriv) %>% pivot_wider(
      names_from = age_m, values_from = paste0("v", var),
      names_prefix = paste0("v", var, ""))
}
