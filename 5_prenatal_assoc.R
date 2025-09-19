#---------------------# 
#### LOAD PACKAGES #### 
#---------------------# 

library(RColorBrewer) 
library(tidyverse) 
library(patchwork) 
library(ggeffects)
library(splines)
library(scales)
library(broom) 
library(lme4)

#------------------------------------------#
#### PRENATAL FACTORS & GROWTH FEATURES #### 
#------------------------------------------# 

load("psme_features.RData") 

psme_est_MF <- psme_est_MF %>% filter(
  !is.na(ReboundBMI) & !is.na(PeakBMI) & 
    !is.na(PeakBMIAge) & !is.na(ReboundBMIAge))

summary(psme_est_MF)

# MISSING: m_age 0, m_ht 20, m_wt 61, GA 0, BWT 0, BLT 2

psme_est_MF <- psme_est_MF %>% filter(
  !is.na(m_ht) & !is.na(m_wt) & !is.na(BLT)) # 823

psme_est_MF %>% select(m_age, m_ht, m_wt, GA, BWT, BLT) %>% summarise_all(
  list(mean = ~round(mean(.,na.rm = T), 1), sd = ~round(sd(.,na.rm = T), 1)))

psme_est_MF <- psme_est_MF %>% mutate_at( c("m_ht", "m_wt", "m_age", "GA", "BLT", "BWT"), .funs = list(z = ~ scale(.))) %>% as.data.frame() 

exp <- c("m_ht_z", "m_wt_z", "m_age_z", "GA_z", "BLT_z", "BWT_z") 
out <- c("PHV", "PWV", "PeakBMI", "ReboundBMI", "PeakBMIAge * 12", "ReboundBMIAge * 12") 

(psme_assoc_res <- expand.grid(out, exp) %>% group_by(Var1) %>% rowwise() %>% 
    summarise(frm = paste0(Var1, " ~ sex + ", Var2)) %>% group_by(model_id = row_number(), frm) %>% 
    do(cbind(tidy( lm(.$frm, data = psme_est_MF)))) %>% mutate( 
      lci = estimate - (1.96 * std.error), uci = estimate + (1.96 * std.error)) %>% 
    filter(term != "sexGirls" & term != "(Intercept)") %>% select(-std.error, -std.error, -statistic) %>% 
    as.data.frame()) 

psme_assoc_res$term <- dplyr::recode(
  psme_assoc_res$term, 
  "m_age_z" = "Maternal age", 
  "m_wt_z" = "Maternal weight", 
  "m_ht_z" = "Maternal height", 
  "GA_z" = "Gestational age", 
  "BWT_z" = "Birth weight", 
  "BLT_z" = "Birth length" 
  ) 

psme_assoc_res$term <- factor(
  psme_assoc_res$term, levels=c(
    "Maternal age", 
    "Maternal weight", 
    "Maternal height", 
    "Gestational age", 
    "Birth weight", 
    "Birth length" )) 

psme_assoc_res$outs <- rep(c(
  "a. Peak height velocity", 
  "b. Peak weight velocity", 
  "c. Peak BMI", 
  "d. Rebound BMI", 
  "e. Age at peak BMI", 
  "f. Age at rebound BMI"), 6) 

rm(psme_est_MF, exp, out) 

# PLOT

assoc_plot <- function(data, outcome_filter, y_lab, plot_title) {
  
  data %>%
    filter(outs == outcome_filter) %>% ggplot(
      aes(x = term, y = estimate, ymin = lci, ymax = uci)) +
    geom_pointrange(aes(color = term)) + geom_hline(
      yintercept = 0, linewidth = 0.2, color = "red") +
    coord_flip() + scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    scale_x_discrete(limits = rev) + labs(
      title = plot_title, y = y_lab, x = NULL) +
    theme_classic() + theme(
      axis.ticks.y = element_blank(),
      legend.position = "none")
}

(p1 <- assoc_plot(psme_assoc_res, "a. Peak height velocity", "cm/yr", "a. Peak height velocity"))
(p2 <- assoc_plot(psme_assoc_res, "b. Peak weight velocity", "kg/yr", "b. Peak weight velocity"))
(p3 <- assoc_plot(psme_assoc_res, "c. Peak BMI", "kg/m²", "c. Peak BMI"))
(p4 <- assoc_plot(psme_assoc_res, "d. Rebound BMI", "kg/m²", "d. Rebound BMI"))
(p5 <- assoc_plot(psme_assoc_res, "e. Age at peak BMI", "months", "e. Age at peak BMI"))
(p6 <- assoc_plot(psme_assoc_res, "f. Age at rebound BMI", "months", "f. Age at rebound BMI"))

graphics.off() 
grDevices::cairo_pdf("res/fig.prenatal_assoc.pdf", height = 3, width = 13.5) 
(p1 | p2 | p3 | p4 | p5 | p6) + plot_layout( 
  guides = "collect", axis_titles = "collect", axes = "collect") & 
  theme( legend.position = 'none') 
dev.off()

# ESTIMATES

psme_assoc_res$res <- paste0(
  round(psme_assoc_res$estimate, 2), " (",
  round(psme_assoc_res$lci, 2), " to ",
  round(psme_assoc_res$uci, 2), ")")

psme_assoc_res %>% write.csv("res/psme_assoc_res.csv")

rm(list = ls())

#-------------------------------------------------#
#### MATERNAL ETHNICITY & HAD/WAD TRAJECTORIES ####
#-------------------------------------------------#

load("psme_features.RData")

psme_est_MF$m_ethn <- dplyr::recode(
  psme_est_MF$m_ethn, 
  "1" = "Chinese", 
  "2" = "Malay", 
  "3" = "Indian"
)

psme_est_MF$m_ethn <- as.factor(psme_est_MF$m_ethn)

table(psme_est_MF$m_ethn)

ethn_HAD_dat <- psme_est_MF %>% select(
  id, sex, m_ethn, HAD_1m, HAD_6m, HAD_12m, HAD_18m, HAD_24m, 
  HAD_30m, HAD_36m, HAD_42m, HAD_48m, HAD_54m, HAD_60m) %>% pivot_longer(
    cols = starts_with("HAD"), names_to = "age", values_to = "had") %>% 
  mutate(t = ifelse(age == "HAD_1m", 1, ifelse(
    age == "HAD_6m", 2, ifelse(
      age == "HAD_12m", 3, ifelse(
        age == "HAD_18m", 4, ifelse(
          age == "HAD_24m", 5, ifelse(
            age == "HAD_30m", 6, ifelse(
              age == "HAD_36m", 7, ifelse(
                age == "HAD_42m", 8, ifelse(
                  age == "HAD_48m", 9, ifelse(
                    age == "HAD_54m", 10, 11
                  )))))))))))

ethn_WAD_dat <- psme_est_MF %>% select(
  id, sex, m_ethn, WAD_1m, WAD_6m, WAD_12m, WAD_18m, WAD_24m, 
  WAD_30m, WAD_36m, WAD_42m, WAD_48m, WAD_54m, WAD_60m) %>% pivot_longer(
    cols = starts_with("WAD"), names_to = "age", values_to = "wad") %>% 
  mutate(t = ifelse(age == "WAD_1m", 1, ifelse(
    age == "WAD_6m", 2, ifelse(
      age == "WAD_12m", 3, ifelse(
        age == "WAD_18m", 4, ifelse(
          age == "WAD_24m", 5, ifelse(
            age == "WAD_30m", 6, ifelse(
              age == "WAD_36m", 7, ifelse(
                age == "WAD_42m", 8, ifelse(
                  age == "WAD_48m", 9, ifelse(
                    age == "WAD_54m", 10, 11
                  )))))))))))

rm(psme_est_MF)

#### ME MODELS

# HAD

HAD_mod <- lmer(
  had ~ (ns(t, df = 3) * m_ethn) + sex + ( 1 | id ), data = ethn_HAD_dat
)

HAD_mod_P <- ggemmeans(HAD_mod, terms = c('m_ethn [all]', 't [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]'))
HAD_mod_P$out <- "a. Height difference, cm"

# WAD

WAD_mod <- lmer(
  wad ~ (ns(t, df = 3) * m_ethn) + sex + ( 1 | id ), data = ethn_WAD_dat
)

WAD_mod_P <- ggemmeans(WAD_mod, terms = c('m_ethn [all]', 't [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]'))
WAD_mod_P$out <- "b. Weight difference, kg"

# PLOT

mods_P <- bind_rows(HAD_mod_P, WAD_mod_P)

mods_P$group <- dplyr::recode(
  mods_P$group, 
  "1" = "1", 
  "1.5" = "3", 
  "2" = "6", 
  "2.5" = "9", 
  "3" = "12", 
  "3.5" = "15", 
  "4" = "18", 
  "4.5" = "21", 
  "5" = "24", 
  "5.5" = "27", 
  "6" = "30", 
  "6.5" = "33", 
  "7" = "36", 
  "7.5" = "39", 
  "8" = "42", 
  "8.5" = "45", 
  "9" = "48",
  "9.5" = "51",
  "10" = "54",
  "10.5" = "57",
  "11" = "60"
  
)

mods_P$group <- factor(
  mods_P$group, levels=c(
    "1",
    "3",
    "6",
    "9",
    "12",
    "15",
    "18",
    "21",
    "24",
    "27",
    "30",
    "33",
    "36",
    "39",
    "42",
    "45",
    "48",
    "51",
    "54",
    "57",
    "60"
  ))

graphics.off()
grDevices::cairo_pdf("res/fig_hwad_ethn.pdf", height = 5, width = 4)

ggplot(mods_P, aes(x = group, y = predicted, group = x, col = x, fill = x)) + theme_classic() + 
  geom_line(linewidth = 1.5) + geom_ribbon(aes(ymin = conf.low, ymax = conf.high), colour = NA, alpha = 0.15) +
  scale_color_brewer(palette = "Set2") +  ylab("Mean difference relative to WHO standards ") + 
  scale_x_discrete(breaks=c(1, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60)) +
  xlab("age, months") + facet_wrap(. ~ out, strip.position = "top", ncol = 1, scales = "free") + theme(
    legend.position = "bottom", legend.title = element_blank(), strip.background = element_blank(),
    strip.text = element_text(face="bold", hjust = 0), axis.ticks.y = element_blank()) + guides(
      override.aes = list(size = 1))

dev.off()

rm(list = ls())
