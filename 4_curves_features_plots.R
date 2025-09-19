
#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(corrplot)
library(rlang)
library(scales)
library(sitar)
library(psme)

#----------------------------------------#
#### PLOT PSME GROWTH/VELOCITY CURVES ####
#----------------------------------------#

load("gusto_psme_preds.RData")
load("gusto_psme_mods.RData")

age.M <- with(dat_M, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))
age.F <- with(dat_F, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))

# GROWTH CURVES

graphics.off()
grDevices::cairo_pdf("res/fig_psme_dist_curves.pdf", height = 5.5, width = 5)
par(mfrow = c(3, 2), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0), mgp = c(2.5, 1, 0))

plot_params <- list(
  list(data = psme_pred_M, pop_data = psme_ht_M, age_data = age.M, y_col = "cht", ylim = c(35, 165), y_label = "cm", title = 'a. Height - boys'),
  list(data = psme_pred_F, pop_data = psme_ht_F, age_data = age.F, y_col = "cht", ylim = c(35, 165), y_label = "cm", title = 'b. Height - girls'),
  list(data = psme_pred_M, pop_data = psme_wt_M, age_data = age.M, y_col = "cwt", ylim = c(0, 80), y_label = "kg", title = 'c. Weight - boys'),
  list(data = psme_pred_F, pop_data = psme_wt_F, age_data = age.F, y_col = "cwt", ylim = c(0, 80), y_label = "kg", title = 'd. Weight - girls'),
  list(data = psme_pred_M, pop_data = psme_bmi_M, age_data = age.M, y_col = "cbmi", ylim = c(8, 35), y_label = "kg/m²", title = 'e. BMI - boys'),
  list(data = psme_pred_F, pop_data = psme_bmi_F, age_data = age.F, y_col = "cbmi", ylim = c(8, 35), y_label = "kg/m²", title = 'f. BMI - girls')
)

for (p in plot_params) {
  mplot(x = age^2, y = p$data[[p$y_col]], id = id, data = p$data, col = id, las = 1,
        ylim = p$ylim, xlab = 'age', ylab = p$y_label)
  lines(p$age_data^2, xTraj(p$pop_data, p$age_data)$pop, type = "l", col = 'black', lwd = 3)
  title(p$title, adj = 0, line = 0.6)
}

dev.off()

# VELOCITY CURVES

calculate_individual_velocity <- function(data, value_col) {
  data %>% select(id, age, !!dplyr::sym(value_col)) %>% mutate(age = age^2) %>% group_by(id) %>% 
    reframe(dY = diff(!!dplyr::sym(value_col)) / diff(age), dX = rowMeans(embed(age, 2))) %>%
    ungroup()
}

psme_htV_M <- calculate_individual_velocity(psme_pred_M, "cht")
psme_wtV_M <- calculate_individual_velocity(psme_pred_M, "cwt")
psme_htV_F <- calculate_individual_velocity(psme_pred_F, "cht")
psme_wtV_F <- calculate_individual_velocity(psme_pred_F, "cwt")

graphics.off()
grDevices::cairo_pdf("res/fig_psme_vel_curves.pdf", height = 5, width = 6)
par(mfrow = c(2, 2), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0), mgp = c(2.5, 1, 0))

plot_params <- list(
  list(data = psme_htV_M, title = 'a. Height velocity - boys', y_label = 'cm/year', ylim = c(0, 135)),
  list(data = psme_htV_F, title = 'b. Height velocity - girls', y_label = 'cm/year', ylim = c(0, 135)),
  list(data = psme_wtV_M, title = 'c. Weight velocity - boys', y_label = 'kg/year', ylim = c(-5, 47)),
  list(data = psme_wtV_F, title = 'd. Weight velocity - girls', y_label = 'kg/year', ylim = c(-5, 47))
)

for (params in plot_params) {
  mplot(
    x = dX, 
    y = dY, 
    id = id, 
    data = params$data, 
    col = id, 
    las = 1, 
    ylim = params$ylim, 
    xlab = 'age', 
    ylab = params$y_label
  )
  title(params$title, adj = 0, line = 0.6)
}

dev.off()

rm(list = ls())

#-----------------------------------------#
#### PLOT PREDICTED VS OBSERVED VALUES ####
#-----------------------------------------#

load("gusto_ps_dat.RData")
load("gusto_psme_preds.RData")

#dat_M %>% group_by(id) %>% filter(n() == 5) %>% ungroup %>% sample_n(1)
#dat_M %>% group_by(id) %>% filter(n() == 10) %>% ungroup %>% sample_n(1)
#dat_M %>% group_by(id) %>% filter(n() == 15) %>% ungroup %>% sample_n(1)
# N=5: 386
# N=10: 460
# N=15: 142

#dat_F %>% group_by(id) %>% filter(n() == 5) %>% ungroup %>% sample_n(1)
#dat_F %>% group_by(id) %>% filter(n() == 10) %>% ungroup %>% sample_n(1)
#dat_F %>% group_by(id) %>% filter(n() == 15) %>% ungroup %>% sample_n(1)
# N=5: 353
# N=10: 339
# N=15: 205

# MALES

psme_pred_M <- psme_pred_M %>% filter(id == 386 | id == 460 | id == 142) %>% 
  rename(pred_ht = cht, pred_wt = cwt, pred_bmi = cbmi) %>% mutate(
    age = age^2)

psme_pred_M <- dat_M %>% filter(
  id == "386" | id == "460" | id == "142") %>% 
  select(id, age, cht, cwt, cbmi) %>% mutate(id = as.integer(id)) %>% 
  full_join(psme_pred_M) %>% mutate(Nobs = ifelse(
    id == 386, "5 points", ifelse(
      id == 460, "10 points", "15 points")
    ), sex = "Boys")

# FEMALES

psme_pred_F <- psme_pred_F %>% filter(id == 352 | id == 339 | id == 205) %>% 
  rename(pred_ht = cht, pred_wt = cwt, pred_bmi = cbmi) %>% mutate(
    age = age^2)

psme_pred_F <- dat_F %>% filter(
  id == "352" | id == "339" | id == "205") %>% 
  select(id, age, cht, cwt, cbmi) %>% mutate(id = as.integer(id)) %>% 
  full_join(psme_pred_F) %>% mutate(Nobs = ifelse(
    id == 352, "5 points", ifelse(
      id == 339, "10 points", "15 points")
      ), sex = "Girls")

# COMBINE & PLOT

psme_pred_MF <- bind_rows(psme_pred_M, psme_pred_F)

psme_pred_MF$Nobs <- factor(
  psme_pred_MF$Nobs, levels=c(
    "5 points",
    "10 points",
    "15 points"
  ))

(ht_plot_MF <- ggplot(
  psme_pred_MF, aes(x = age, col = sex)) + theme_classic() + 
    geom_line(aes(y = pred_ht)) + geom_point(aes(y = cht), size = 0.5) +
    facet_wrap(. ~ Nobs, strip.position = "top", ncol = 4) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle("a. Height") + ylab("cm") + theme(
      legend.position = "bottom", strip.background = element_blank(),
      legend.title = element_blank(), strip.text = element_text(face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=3))))

(wt_plot_MF <- ggplot(
  psme_pred_MF, aes(x = age, col = sex)) + theme_classic() + 
    geom_line(aes(y = pred_wt)) + geom_point(aes(y = cwt), size = 0.5) +
    facet_wrap(. ~ Nobs, strip.position = "top", ncol = 4) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle("b. Weight") + ylab("kg") + theme(
      legend.position = "bottom", strip.background = element_blank(),
      legend.title = element_blank(), strip.text = element_text(face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=3))))

(bmi_plot_MF <- ggplot(
  psme_pred_MF, aes(x = age, col = sex)) + theme_classic() + 
    geom_line(aes(y = pred_bmi)) + geom_point(aes(y = cbmi), size = 0.5) +
    facet_wrap(. ~ Nobs, strip.position = "top", ncol = 4) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10)) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle("c. BMI") + ylab("kg/m²") + theme(
      legend.position = "bottom", strip.background = element_blank(),
      legend.title = element_blank(), strip.text = element_text(face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=3))))

graphics.off()
grDevices::cairo_pdf("res/fig_pred_obs_MF.pdf", height = 5.8, width = 5)
(ht_plot_MF / wt_plot_MF / bmi_plot_MF) + plot_layout(
  guides = "collect", axis_titles = "collect", axes = "collect") & theme(
    legend.position = 'bottom')
dev.off()

rm(list = ls())

#---------------------------------#
#### BOXPLOTS GROWTH FEATURES  ####
#--------====---------------------#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% select(
  sex, PWV, PHV, PeakBMI, PeakBMIAge, ReboundBMI, ReboundBMIAge)

psme_est_MF <- psme_est_MF %>% 
  mutate(PeakBMIAge = PeakBMIAge * 12) %>% pivot_longer(cols = c(
    PWV, PHV, PeakBMI, PeakBMIAge, ReboundBMI, ReboundBMIAge), 
    names_to = c("outs"), values_to = "count")

psme_est_MF$outs <- dplyr::recode(
  psme_est_MF$outs, 
  PHV = "a. Peak height velocity (cm/yr)",
  PWV = "b. Peak weight velocity (kg/yr)",
  PeakBMI = "c. Peak BMI (kg/m²)",
  ReboundBMI = "d. Rebound BMI (kg/m²)",
  PeakBMIAge = "e. Age at peak BMI (months)",
  ReboundBMIAge = "f. Age at rebound BMI (years)"
)

graphics.off()
grDevices::cairo_pdf("res/fig_boxplot_features_raw.pdf", height = 5.5, width = 5.5)

ggplot(psme_est_MF) + theme_classic() + 
  geom_boxplot(aes(x = count, col = sex), outlier.size = 0.2) + 
  scale_color_brewer(palette = "Set1") + 
  facet_wrap(outs ~ ., scales = "free", ncol = 2) + coord_flip() + theme(
    legend.title = element_blank(), legend.position = "right", 
    axis.title = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), strip.background = element_blank(),
    strip.text = element_text(face="bold", hjust = 0)) + guides(
      colour = guide_legend(override.aes = list(
        size = 5)))

dev.off()

psme_est_MF %>% group_by(outs, sex) %>% summarise_all(list(
  median = ~round(median(.,na.rm = T), 1),
  q1 = ~round(quantile(., na.rm = T, probs = 0.25), 1),
  q3 = ~round(quantile(., na.rm = T, probs = 0.75), 1),
  mean = ~round(mean(.,na.rm = T), 1),
  sd = ~round(sd(.,na.rm = T), 1))) %>% 
  write.csv("res/features.csv")

rm(list = ls())

#----------------------------------------#
#### GROWTH FEATURES CORRELATION PLOT ####
#----------------------------------------#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% select(
  sex, PWV, PHV, PeakBMI, PeakBMIAge, ReboundBMI, ReboundBMIAge) %>% drop_na()

psme_est_MF <- psme_est_MF %>% rename(
  "Peak weight velocity" = PWV,
  "Peak height velocity" = PHV,
  "Peak BMI" = PeakBMI,
  "Rebound BMI" = ReboundBMI,
  "Age at peak BMI" = PeakBMIAge,
  "Age at rebound BMI" = ReboundBMIAge
)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

graphics.off()
grDevices::cairo_pdf("res/fig_corplot_features.pdf", height = 7.5, width = 7.5)
par(mfrow = c(2,1), mar = c(4, 3.6, 1.6, 0.5), oma = c(0, 4, 0, 0))

corrplot(cor(psme_est_MF %>% filter(sex == "Boys") %>% select(-sex)), 
         method = "color", col=col(200),  
         type = "upper", order = "alphabet",
         tl.col = "black", tl.srt = 45,
         diag=FALSE,
         addCoef.col = "black"
)

mtext("a. Boys", side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)

corrplot(cor(psme_est_MF %>% filter(sex == "Girls") %>% select(-sex)), 
         method = "color", col=col(200),  
         type = "upper", order = "alphabet",
         tl.col = "black", tl.srt = 45,
         diag=FALSE,
         addCoef.col = "black"
)

mtext("b. Girls", side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)

dev.off()

rm(list = ls())