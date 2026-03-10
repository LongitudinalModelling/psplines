
#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(RColorBrewer)
library(ggcorrplot)
library(ggeffects)
library(tidyverse)
library(patchwork)
library(corrplot)
library(splines)
library(scales)
library(lme4)

#---------------------------------#
#### BOXPLOTS GROWTH FEATURES  ####
#--------====---------------------#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% select(
  sex, PWV, PHV, PeakBMI, PeakBMIAge, ReboundBMI, ReboundBMIAge)

psme_est_MF <- psme_est_MF %>% pivot_longer(cols = c(
    PWV, PHV, PeakBMI, PeakBMIAge, ReboundBMI, ReboundBMIAge), 
    names_to = c("outs"), values_to = "count")

psme_est_MF$outs <- dplyr::recode(
  psme_est_MF$outs, 
  PHV = "a. Peak height velocity (cm/mo)",
  PWV = "b. Peak weight velocity (g/mo)",
  PeakBMI = "c. Peak BMI (kg/m²)",
  ReboundBMI = "d. Rebound BMI (kg/m²)",
  PeakBMIAge = "e. Age at peak BMI (mo)",
  ReboundBMIAge = "f. Age at rebound BMI (y)"
)

graphics.off()
grDevices::cairo_pdf("res/fig_boxplot_features_raw.pdf", height = 5.5, width = 5)

ggplot(psme_est_MF) + theme_classic() + 
  geom_boxplot(aes(x = count, col = sex), outlier.size = 0.2) + 
  scale_color_brewer(palette = "Set1") + 
  facet_wrap(outs ~ ., scales = "free", ncol = 2) + coord_flip() + theme(
    legend.title = element_blank(), legend.position = "bottom", 
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

#--------------------------------------#
#### BOXPLOTS/CORR INFANT VELOCITY  ####
#--------====--------------------------#

#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% 
  select(sex, vcwt1, vcwt6, vcwt12, vcwt24, vcht1, vcht6, vcht12, vcht24) %>% 
  pivot_longer(
    cols = c(vcwt1, vcwt6, vcwt12, vcwt24, vcht1, vcht6, vcht12, vcht24),
    names_to = "outs",
    values_to = "count") %>% mutate(
    domain = case_when(
      grepl("vcht", outs) ~ "a. Height velocity (cm/mo)",
      grepl("vcwt", outs) ~ "b. Weight velocity (g/mo)"
    ),
    time = recode(
      outs,
      vcht1 = "1m", vcht6 = "6m", vcht12 = "12m", vcht24 = "24m",
      vcwt1 = "1m", vcwt6 = "6m", vcwt12 = "12m", vcwt24 = "24m"
    )
  )

psme_est_MF$time <- factor(
  psme_est_MF$time, levels=c(
    "1m",
    "6m",
    "12m",
    "24m"
  ))

(p_box <- ggplot(psme_est_MF, aes(x = time, y = count, col = sex)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~ domain, scales = "free_y") +
  theme_classic() + theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", , size = 12, hjust = 0)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5))))

psme_est_MF %>% group_by(outs, sex) %>% select(count) %>% summarise_all(list(
  median = ~round(median(.,na.rm = T), 1),
  q1 = ~round(quantile(., na.rm = T, probs = 0.25), 1),
  q3 = ~round(quantile(., na.rm = T, probs = 0.75), 1),
  mean = ~round(mean(.,na.rm = T), 1),
  sd = ~round(sd(.,na.rm = T), 1))) %>% 
  write.csv("res/vel_features.csv")

#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% 
  select(sex, vcwt1, vcwt6, vcwt12, vcwt24, vcht1, vcht6, vcht12, vcht24)

psme_est_MF <- psme_est_MF %>% rename(
  "Height: 1m" = vcht1,
  "Height: 6m" = vcht6,
  "Height: 12m" = vcht12,
  "Height: 24m" = vcht24,
  "Weight: 1m" = vcwt1,
  "Weight: 6m" = vcwt6,
  "Weight: 12m" = vcwt12,
  "Weight: 24m" = vcwt24
)

df <- psme_est_MF %>%
  group_by(sex) %>%
  nest() %>%
  mutate(
    plot_title = case_when(
      sex == "Boys" ~ "c. Growth velocity correlations: boys",
      sex == "Girls" ~ "d. Growth velocity correlations: girls",
      TRUE ~ sex # Fallback for any other 'sex' value
    ),
    plot = map2(data, plot_title, ~ggcorrplot(
      cor(select_if(.x, is.numeric)), 
      type = "lower", 
      lab = TRUE, 
      show.legend = FALSE
    ) + ggtitle(.y) +
      theme(plot.title = element_text(face = "bold", hjust = 2)))
  )

(p_cor <- df$plot %>% wrap_plots(ncol = 2))

#### combine and save

p_box2 <- (p_box + plot_spacer()) + plot_layout(widths = c(2, 0.05))

graphics.off()
grDevices::cairo_pdf("res/fig_boxplot_features_corr.pdf", height = 7, width = 9)

(p_box2 / p_cor) + plot_layout(heights = c(1, 1.25))

dev.off()

rm(list = ls())
