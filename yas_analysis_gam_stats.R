library(tidyverse)
library(lmerTest)
library(ggpubr)
library(corrplot)
library(ez)
library(Hmisc)
library(car)
library(emmeans)
source("~/Documents/R functions/corstars.R")

# load data ####
df_gam <- read.table(file = "~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/gam/df_gam_itemLevel_dependentSE.txt")

df_ppt <- read_delim("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/ppt_demos.txt", delim = "\t") %>%
  filter(exclude != 1) %>%
  select(ppt, CHI, BILINGUAL, ENG, 
         chinese_start_age, chinese_daily_use, chinese_speaking, chinese_listening, chinese_reading, chinese_writing, chinese_overall, chinese_listeningspeaking,
         english_start_age, english_daily_use, english_speaking, english_listening, english_reading, english_writing, english_overall, english_listeningspeaking,
         music_formal, music_informal, music_total) %>%
  droplevels()


# add group info and pick participants
df_gam_complete <- df_gam %>%
  inner_join(df_ppt, by = "ppt") %>%
  mutate(group = case_when(
    ENG==1 ~ "ENG", # ENG mean English speakers excluding SEA and CBC
    CHI=="1" ~ "CHI",
    BILINGUAL=="1" ~ "BILINGUAL",
    .default = "ELSE"
  )) %>%
  # select group and conditions
  filter(group %in% c("ENG", "CHI", "BILINGUAL")) %>%
  filter(condition %in% c("ChEn-devi", "EnCh-devi")) %>%
  droplevels() %>%
  # convert hasPeak to 1 and 0
  mutate(hasPeak_numeric = as.numeric(hasPeak))


# plotting hasPeak
ggplot(data = df_gam_complete) +
  geom_bar(aes(x = group, y = hasPeak_numeric, fill = condition),
           stat = "summary", fun = mean, position=position_dodge()) +
  coord_cartesian(ylim = c(0.8, 1)) +
  theme_bw()

no_peak_ppts <- df_gam_complete %>%
  group_by(ppt) %>%
  dplyr::summarize(peak_mean = mean(hasPeak_numeric)) %>%
  ungroup() %>%
  filter(peak_mean < 1) %>%
  droplevels()

df_gam_hasPeak <- df_gam_complete %>%
  # filter(hasPeak == "TRUE") %>%
  filter(hasPeak == "TRUE" & !(ppt %in% no_peak_ppts$ppt) ) %>%
  # filter(half_area_latency >= 150 & half_area_latency <= 600) %>%
  droplevels()

# check sample size
df_gam_hasPeak %>% group_by(group, condition) %>% dplyr::summarize(n = n())

# plotting and stats ####

## plot data ####
df_plot <- df_gam_hasPeak %>%
  # filter(condition %in% c("ChEn-devi")) %>%
  filter(condition %in% c("ChEn-devi", "EnCh-devi")) %>%
  mutate(
    condition = factor(condition, levels = c("ChEn-devi", "EnCh-devi")),
    group = factor(group, levels = c("ENG", "BILINGUAL", "CHI"))
  ) %>%
  droplevels()

## stats data ####
df_mod <- df_gam_hasPeak %>%
  # filter(condition %in% c("ChEn-devi")) %>%
  filter(condition %in% c("ChEn-devi", "EnCh-devi")) %>%
  mutate(
    condition = factor(condition, levels = c("ChEn-devi", "EnCh-devi")),
    group = factor(group, levels = c("ENG", "BILINGUAL", "CHI"))
         ) %>%
  droplevels()

# set contrast
en_ch <- c(-1/2, 1/2)
contrasts(df_mod$condition) <- cbind(en_ch)

LANG_AUDI <- c(-1/3, -1/3, 2/3)
MONO_BILI <- c(-1/2, 1/2, 0)
contrasts(df_mod$group) <- cbind(LANG_AUDI, MONO_BILI)

## traditional erp ####
fig <-
  ggplot(data = df_plot) +
  geom_boxplot(aes(x = group, y = trad_erp, fill = condition)) +
  theme_bw()
print(fig)
# ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_trad_erp.png")

# mod <- lm(trad_erp ~ group, data = df_mod)
mod <- lmer(trad_erp ~ group*condition + (1|ppt), data = df_mod)
summary(mod)
Anova(mod, type="III")
pair_comp <- emmeans(mod, pairwise ~ group | condition)
test(pair_comp)



## gam erp ####
fig <-
  ggplot(data = df_plot) +
  geom_boxplot(aes(x = group, y = gam_erp, fill = condition)) +
  theme_bw()
print(fig)
# ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_trad_erp.png")

# mod <- lm(gam_erp ~ group, data = df_mod)
mod <- lmer(gam_erp ~ group*condition + (1|ppt), data = df_mod)
summary(mod)
Anova(mod, type = "III")
pair_comp <- emmeans(mod, pairwise ~ group | condition)
test(pair_comp)


## modeled area ####
fig <- 
  ggplot(data = df_plot) +
  geom_boxplot(aes(x = group, y = area, fill = condition)) +
  theme_bw()
print(fig)

# mod <- lm(area ~ group, data = df_mod)
mod <- lmer(area ~ group*condition + (1|ppt), data = df_mod)
summary(mod)
Anova(mod, type="III")
# pair_comp <- emmeans(mod, pairwise ~ group)
pair_comp <- emmeans(mod, pairwise ~ group | condition)
test(pair_comp)

## modeled peak ####
fig <- 
  ggplot(data = df_plot) +
  geom_boxplot(aes(x = group, y = peak_height, fill = condition)) +
  theme_bw()
print(fig)

# mod <- lm(peak_height ~ group, data = df_mod)
mod <- lmer(peak_height ~ group*condition + (1|ppt), data = df_mod)
summary(mod)
Anova(mod)
# pair_comp <- emmeans(mod, pairwise ~ group)
pair_comp <- emmeans(mod, pairwise ~ group | condition)
test(pair_comp)
pair_comp <- emmeans(mod, pairwise ~ group)
test(pair_comp)

## normalized modeled peak ####
fig <- 
  ggplot(data = df_plot) +
  geom_boxplot(aes(x = group, y = NMP, fill = condition)) +
  theme_bw()
print(fig)
# mod <- lm(NMP ~ group, data = df_mod)
mod <- lmer(NMP ~ group*condition + (1|ppt), data = df_mod)
summary(mod)
Anova(mod, type="III")
# pair_comp <- emmeans(mod, pairwise ~ group)
pair_comp <- emmeans(mod, pairwise ~ group | condition)
test(pair_comp)

## modeled fractional area latency ####
fig <- 
  ggplot(data = df_plot) +
  geom_boxplot(aes(x = group, y = half_area_latency, fill = condition)) +
  theme_bw()
print(fig)
# mod <- lm(half_area_latency ~ group, data = df_mod)
mod <- lmer(half_area_latency ~ group*condition + (1|ppt), data = df_mod)
summary(mod)
Anova(mod, type="III")
# pair_comp <- emmeans(mod, pairwise ~ group)
pair_comp <- emmeans(mod, pairwise ~ group | condition)
test(pair_comp)

## modeled peak latency ####
fig <- 
  ggplot(data = df_plot) +
  geom_boxplot(aes(x = group, y = peak_time, fill = condition)) +
  theme_bw()
print(fig)
# mod <- lm(peak_time ~ group, data = df_mod)
mod <- lmer(peak_time ~ group*condition + (1|ppt), data = df_mod)
summary(mod)
Anova(mod, type="III")
# pair_comp <- emmeans(mod, pairwise ~ group)
pair_comp <- emmeans(mod, pairwise ~ group | condition)
test(pair_comp)



# Brain-behavioral correlation ####
# read in reading data
df_rating_raw <- read_delim("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/ratings.txt")
df_rating <- na.omit(df_rating_raw) %>%
  inner_join(df_gam_hasPeak, by = "ppt") %>%
  mutate(language = case_when(
    group %in% c("ENG", "CHI") ~ "mono",
    group %in% c("BILINGUAL") ~ "bili"
  )) %>%
  mutate(
    group = factor(group, levels = c("ENG", "CHI", "BILINGUAL")),
    condition = factor(condition, levels = c("ChEn-devi", "EnCh-devi")),
    language = factor(language, levels = c("mono", "bili"))) #%>% filter(half_area_latency>=150 & half_area_latency<=600)


df_cor <- df_rating %>%
  filter(group %in% c("CHI") & condition == "EnCh-devi") %>%
  # filter(condition == "ChEn-devi") %>%
  # filter(music_formal>0 & music_formal<25) %>%
  droplevels() %>%
  select(
    c(
      "familiarity_chinese",
      "likelihood_chinese",
      # "familiarity_indian",
      # "likelihood_indian",
      "trad_erp",
      "gam_erp",
      "area",
      "peak_height",
      "NMP",
      "peak_time",
      "half_area_latency",
      # for BILINGUAL group:
      # "chinese_start_age",
      # "chinese_daily_use",
      # "chinese_speaking",
      # "chinese_listening",
      # "chinese_reading",
      # "chinese_writing",
      # "chinese_overall",
      # "chinese_listeningspeaking",
      # # for CHI group
      # "english_start_age",
      # "english_daily_use",
      # "english_speaking",
      # "english_listening",
      # "english_reading",
      # "english_writing",
      "english_overall",
      "english_listeningspeaking"
      # # music
      # "music_formal",
      # "music_informal",
      # "music_total"
    ))
corstars(df_cor)


# cor(df_cor, use = "complete.obs")
cor_matrix <- cor(df_cor, use = "pairwise.complete.obs")


# matrix of the p-value of the correlation
p.mat <- cor.mtest(df_cor)

# Extract upper triangle
p_vals <- p.mat[upper.tri(p.mat)]

# Adjust for multiple testing
p_vals_adj <- p.adjust(p_vals, method = "holm")

p.mat[upper.tri(p.mat)] <- p_vals_adj

# plotting
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

# significance marker
corrplot(cor_matrix, method="color", col = col(200),
         type="full", 
         # order="hclust",
         tl.cex = 0.7, tl.col="black", tl.srt = 45, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig",
         pch.col = "grey20", pch.cex = 0.9,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE
)

# correlation coefficient
corrplot(cor_matrix, method="color", col = col(200),
         addCoef.col = "black", number.cex = 0.7,
         type="full", 
         # order="hclust",
         tl.cex = 0.7, tl.col="black", tl.srt = 45, # Text label color and rotation
         # Combine with significance
         # p.mat = p.mat.adj, sig.level = -1, insig = "p-value",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE
)


# correlation plot ####
df_cor_plot <- df_rating %>%
  pivot_longer(cols = c(
    # likelihood_chinese, 
    # familiarity_chinese,
    chinese_overall,
    english_overall,
    ), names_to = "exposure_type", values_to = "rating") %>%
  filter(condition %in% c("ChEn-devi", "EnCh-devi")) %>%
  # filter(half_area_latency>=150 & half_area_latency<=600) %>%
  # filter(group %in% c("ENG", "BILINGUAL")) %>%
  # filter(NMP < -1) %>%
  droplevels()


df_cor_plot <- df_rating %>%
  filter(music_formal>0) %>%
  # filter(group %in% c("CHI")) %>%
  pivot_longer(cols = c(music_formal), names_to = "exposure_type", values_to = "rating") %>%
  filter(condition %in% c("ChEn-devi", "EnCh-devi")) %>%
  mutate(rating = as.numeric(rating)) %>%
  droplevels()


fig <-
  ggplot(data = df_cor_plot,
         mapping = aes(x = rating, y = trad_erp, 
                       group = group, color = group, fill = group
                       )) +
  facet_grid(condition ~ exposure_type) +
  geom_point(size = 1, alpha = 0.8) +
  # coord_cartesian(ylim = c(-7, 7),
  #                 xlim = c(0, 10)) +
  # scale_x_continuous(breaks = seq(1,10, 1)) +
  geom_smooth(method = "lm",
              se = TRUE,
              linewidth = 0.5,
              alpha = 0.1,
              na.rm = TRUE) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
           label.x.npc = 0, #adjust the label in x axis
           label.y.npc = 0.2, #adjust the label in y axis
           size = 3) +
  theme_bw() +
  theme(legend.title = element_blank())
print(fig)
# save 
ggsave(plot = fig, width = 12, height = 7, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/correlation/NMP.png")

# regression

# set contrast
levels(df_rating$condition)
en_ch <- c(-1/2, 1/2)
contrasts(df_rating$condition) <- cbind(en_ch)

levels(df_rating$group)
MONO_BILI <- c(-1/3, -1/3, 2/3)
ENG_CHI <- c(-1/2, 1/2, 0)
contrasts(df_rating$group) <- cbind(MONO_BILI, ENG_CHI)

mod <- lmer(NMP ~ group * condition * scale(likelihood_chinese) + (1|ppt), data = df_rating)
summary(mod)
Anova(mod, type = "III")

pair_comp <- emtrends(mod, pairwise ~ group | condition, var = "likelihood_chinese")
# slope significance
test(pair_comp)
# slope difference
pair_comp
# effect size
eff_size(pair_comp$emtrends, sigma = sigma(mod), edf = df.residual(mod), method = "identity")
