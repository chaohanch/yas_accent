library(tidyverse)
library(lmerTest)
library(ggpubr)
library(corrplot)
library(ez)
library(Hmisc)
source("~/Documents/R functions/corstars.R")

df_gam <- read.table(file = "~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/gam/df_gam.txt")

df_ppt <- read_delim("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/ppt_demos.txt", delim = "\t") %>%
  filter(exclude != 1) %>%
  select(ppt, grouping1, english_start_age, english_daily_use, english_speaking, english_listening, english_reading, english_writing) %>%
  droplevels()


# add group info and pick participants
df_gam_complete <- df_gam %>%
  inner_join(df_ppt, by = "ppt") %>%
  mutate(group = case_when(
    grouping1=="EN" ~ "ENG", # ENG mean English speakers excluding SEA and CBC
    grouping1=="CH" ~ "CHI",
    grouping1=="SE" ~ "SEA",
    grouping1=="CBC" ~ "CBC",
    .default = "ELSE"
  )) %>%
  # select group and conditions
  filter(group %in% c("ENG", "CHI", "SEA", "CBC")) %>%
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
  filter(hasPeak == "TRUE" & !(ppt %in% no_peak_ppts$ppt) ) %>%
  droplevels()

df_mod <- df_gam_hasPeak %>%
  filter(condition %in% c("ChEn-devi", "EnCh-devi")) %>%
  droplevels() %>%
  mutate(condition = factor(condition, levels = c("EnCh-devi", "ChEn-devi")),
         group = factor(group, levels = c("CHI", "ENG")))
ch_en <- c(-1/2, 1/2)
CHIgroup_ENGgroup<- c(-1/2, 1/2)
# CHgroup_Natives<- c(-2/3, 1/3, 1/3)
# ENgroup_CBCgroup <- c(0, -1/2, 1/2)
contrasts(df_mod$condition) <- cbind(ch_en)
contrasts(df_mod$group) <- cbind(CHIgroup_ENGgroup)

## plot traditional erps and stats ####
fig <-
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = trad_erp1, fill = condition)) +
  theme_bw()
print(fig)
ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_trad_erp1.png")

# stats
summary(mod <- lmer(trad_erp1 ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
        dv = trad_erp1,
        wid = ppt,
        within = condition,
        between = group,
        type = 3,
        detailed = TRUE)

fig <-
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = trad_erp2, fill = condition)) +
  theme_bw()
print(fig)
ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_trad_erp2.png")

summary(mod <- lmer(trad_erp2 ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
        dv = trad_erp2,
        wid = ppt,
        within = condition,
        between = group,
        type = 3,
        detailed = TRUE)

## plot gam erps and stats ####
fig <- 
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = gam_erp1, fill = condition)) +
  theme_bw()
print(fig)
ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_gam_erp1.png")

summary(mod <- lmer(gam_erp1 ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
        dv = gam_erp1,
        wid = ppt,
        within = condition,
        between = group,
        type = 3,
        detailed = TRUE)

fig <- 
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = gam_erp2, fill = condition)) +
  theme_bw()
print(fig)
ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_gam_erp2.png")

summary(mod <- lmer(gam_erp2 ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
        dv = trad_erp2,
        wid = ppt,
        within = condition,
        between = group,
        type = 3,
        detailed = TRUE)


## plotting modeled area and stats ####
fig <- 
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = area, fill = condition)) +
  theme_bw()
print(fig)
ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_area.png")

summary(mod <- lmer(area ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
  dv = area,
  wid = ppt,
  within = condition,
  between = group,
  type = 3,
  detailed = TRUE)

## plotting modeled peak amplitude and stats ####
fig <- 
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = peak_height, fill = condition)) +
  theme_bw()
print(fig)
ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_peak_height.png")

summary(lmer(peak_height ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
        dv = peak_height,
        wid = ppt,
        within = condition,
        between = group,
        type = 3,
        detailed = TRUE)

## plotting normalized modeled peak amplitude and stats ####
fig <- 
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = NMP, fill = condition)) +
  theme_bw()
print(fig)
ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_NMP.png")

summary(lmer(NMP ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
        dv = NMP,
        wid = ppt,
        within = condition,
        between = group,
        type = 3,
        detailed = TRUE)

## modeled fractional area latency ####
fig <- 
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = half_area_latency, fill = condition)) +
  theme_bw()
print(fig)
ggsave(plot = fig, width = 5, height = 5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/gam_half_area_latency.png")

summary(lmer(half_area_latency ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
        dv = half_area_latency,
        wid = ppt,
        within = condition,
        between = group,
        type = 3,
        detailed = TRUE)


## modeled peak latency ####
ggplot(data = df_gam_hasPeak) +
  geom_boxplot(aes(x = group, y = peak_time, fill = condition)) +
  theme_bw()
summary(lmer(peak_time ~ condition * group + (1|ppt), data = df_mod))
ezANOVA(data = df_mod,
        dv = peak_time,
        wid = ppt,
        within = condition,
        between = group,
        type = 3,
        detailed = TRUE)


# Brain-behavioral correlation ####
# read in reading data
df_rating_raw <- read_delim("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/ratings.txt")
df_rating <- na.omit(df_rating_raw) %>%
  select(-ppt_group) %>%
  inner_join(df_gam_hasPeak, by = "ppt") %>%
  mutate(language = case_when(
    group %in% c("ENG", "SEA", "CBC") ~ "english_speaker", # ENG mean English speakers excluding SEA and CBC
    group %in% c("CHI") ~ "chinese_speaker",
    .default = "ELSE")) %>%
  mutate(group = factor(group, levels = c("ENG", "SEA", "CBC", "CHI")),
         condition = factor(condition, levels = c("ChEn-devi", "EnCh-devi")),
         language = factor(language, levels = c("english_speaker", "chinese_speaker")))


df_cor <- df_rating %>%
  # filter(group=="ENG" & condition == "ChEn-devi") %>%
  filter(condition == "ChEn-devi") %>%
  droplevels() %>%
  select(
    c(
      "familiarity_chinese",
      "likelihood_chinese",
      "familiarity_indian",
      "likelihood_indian",
      "area",
      "peak_height",
      "NMP",
      "peak_time",
      "half_area_latency",
      # # for CHI speakers:
      # "english_start_age",
      # "english_daily_use",
      # "english_speaking",
      # "english_listening",
      # "english_reading",
      # "english_writing"
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


df_cor_plot <- df_rating %>%
  pivot_longer(cols = familiarity_chinese:likelihood_indian, names_to = "exposure_type", values_to = "rating") %>%
  # filter(group != "CBC") %>%
  ungroup()


fig <-
  ggplot(data = df_cor_plot,
         mapping = aes(x = rating, y = peak_height, group = group, color = group, fill = group)) +
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
levels(df_rating$group)
Native_Nonnative <- c(-1/4, -1/4, -1/4, 3/4)
CHaccent_Else <- c(-1/3, -1/3, 2/3, 0)
INaccent_Else <- c(-1/2, 1/2, 0, 0)
contrasts(df_rating$group) <- cbind(Native_Nonnative, CHaccent_Else, INaccent_Else)

levels(df_rating$condition)
ChDevi_EnDevi <- c(-1/2, 1/2)
contrasts(df_rating$condition) <- cbind(ChDevi_EnDevi)

levels(df_rating$language)
English_Chinese <- c(-1/2, 1/2)
contrasts(df_rating$language) <- cbind(English_Chinese)

mod <- lmer(scale(trad_erp1) ~ group * condition * scale(likelihood_chinese) + (1|condition),
          data = df_rating)
Anova(mod, type = "III")

library(emmeans)
pair_comp <- emtrends(mod, pairwise ~ condition, var = "likelihood_chinese")
# slope significance
test(pair_comp)
# slope difference
pair_comp
# effect size
eff_size(pair_comp$emtrends, sigma = sigma(mod), edf = df.residual(mod), method = "identity")


df_subset <- df_rating[df_rating$language=="english_speaker", ]
mod <- lmer(scale(NMP) ~ condition * (scale(likelihood_chinese) + scale(familiarity_chinese)) + (1|condition),
            data = df_subset)
Anova(mod, type = "III")

df_subset <- df_rating[df_rating$language=="chinese_speaker", ]
mod <- lmer(scale(NMP) ~ condition * (scale(likelihood_chinese) + scale(english_start_age)) + (1|condition),
            data = df_subset)
Anova(mod, type = "III")

mod <- lmer(scale(NMP) ~ condition * scale(likelihood_chinese) + (1|condition),
            data = df_subset)
Anova(mod, type = "III")

library(emmeans)
pair_comp <- emtrends(mod, pairwise ~ condition, var = "likelihood_chinese")
# slope significance
test(pair_comp)
# slope difference
pair_comp
# effect size
eff_size(pair_comp$emtrends, sigma = sigma(mod), edf = df.residual(mod), method = "identity")