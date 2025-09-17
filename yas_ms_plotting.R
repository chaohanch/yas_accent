library(tidyverse)
library(cowplot)
library(ggpubr)
library(magick)

source("Documents/GitHub/yas_accent/yas_ms_theme_function.R")

# demographic ####
df_demo <- read_delim("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/ppt_demos.txt", delim = "\t") %>%
  mutate(group = case_when(
    ENG==1 ~ "ENG", # ENG mean English speakers excluding SEA and CBC
    CHI=="1" ~ "CHI",
    BILINGUAL=="1" ~ "BILINGUAL",
    .default = "ELSE"
  )) %>%
  # select group and conditions
  filter(group %in% c("ENG", "CHI", "BILINGUAL"))

df_demo %>% group_by(Gender) %>% summarize(n = n())
df_demo %>% group_by(group) %>% summarize(n = n())

summary(as.numeric(df_demo$Age...11))
# erp descriptive plot ####

df_erp <- read_delim("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/erp_all.txt", delim = "\t") %>%
  rowwise() %>%
  mutate(amp = mean(Cz, Fz, FC1, FC2, CP1, CP2, C3, C4, Pz)) %>%
  ungroup() %>%
  select(participant, group, condition, time, amp) %>%
  filter(condition %in% c("ChEn-devi", "ChEn-stan", "EnCh-devi", "EnCh-stan")) %>%
  mutate(accent = case_when(
    condition %in% c("ChEn-devi", "EnCh-stan") ~ "SE-accented English",
    condition %in% c("EnCh-devi", "ChEn-stan") ~ "MC-accented English",
    .default = NA)) %>%
  mutate(stim_role = case_when(
    condition %in% c("ChEn-stan", "EnCh-stan") ~ "standard",
    condition %in% c("EnCh-devi", "ChEn-devi") ~ "deviant",
    .default = NA)) %>%
  mutate(
    group = factor(group, levels = c("ENG", "CHI", "BILINGUAL"), labels = c("English speakers", "Mandarin-dominant bilinguals", "English-dominant bilinguals")),
    stim_role = factor(stim_role, levels = c("standard", "deviant")),
    accent = factor(accent, levels = c("SE-accented English", "MC-accented English"))
    )

ylimit <- c(-4, 4)

df_plot <- df_erp %>%
  filter(group == "English speakers") %>%
  droplevels()
fig_english <-
  myfunc_erp_plot(df_plot, "English speakers")
fig_english

df_plot <- df_erp %>%
  filter(group == "English-dominant bilinguals") %>%
  droplevels()
fig_engbiling <-
  myfunc_erp_plot(df_plot, "English-dominant bilinguals")
fig_engbiling

df_plot <- df_erp %>%
  filter(group == "Mandarin-dominant bilinguals") %>%
  droplevels()
fig_manbiling <-
  myfunc_erp_plot(df_plot, "Mandarin-dominant bilinguals")
fig_manbiling

fig <-
  ggdraw() +
  draw_plot(fig_english, x = 0, y = 0, width = 1/3, height = 1) +
  draw_plot(fig_engbiling, x = 1/3, y = 0, width = 1/3, height = 1) +
  draw_plot(fig_manbiling, x = 2/3, y = 0, width = 1/3, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), x = c(0, 1/3, 2/3), y = c(1, 1, 1))
fig

ggsave(plot = fig, width = 9, height = 5.5, units = "in", dpi = 300, filename = "OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/erp_wave.jpg")

# erp permutation plot ####

## English speakers, SE accent ####

sig_times <- df_erp$time[df_erp$time >= 452 & df_erp$time <= 628]
chans <- c('Fz','F3','FC5','FC1','C3','CP5','CP1','Pz','P3','O2','P4','P8','CP6','CP2','Cz','C4','FC2','F4')

### absolute waves ####
df_abs <- read_delim("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/erp_all.txt", delim = "\t") %>%
  mutate( amp = rowMeans(across(all_of(chans))) ) %>%
  select(participant, group, condition, time, amp) %>%
  filter(condition %in% c("ChEn-devi", "EnCh-stan")) %>%
  mutate(stim_role = case_when(
    condition %in% c("ChEn-stan", "EnCh-stan") ~ "standard",
    condition %in% c("EnCh-devi", "ChEn-devi") ~ "deviant",
    .default = NA)) %>%
  mutate(
    group = factor(group, levels = c("ENG", "CHI", "BILINGUAL"), labels = c("English speakers", "Mandarin-dominant bilinguals", "English-dominant bilinguals")),
    stim_role = factor(stim_role, levels = c("standard", "deviant"))
  ) %>%
  filter(group=="English speakers") %>%
  droplevels()

# plot
ylimit <- c(-2, 2)
fig_abs <- myfunc_permute_abs_plot(df_abs, "English speakers", "SE accent")
fig_abs

### difference waves ####
df_diff <- df_abs %>%
  select(-stim_role) %>%
  pivot_wider(names_from = condition, values_from = amp) %>%
  mutate(amp = `ChEn-devi` - `EnCh-stan`,
         condition = "deviant - standard")
  
# mutate(condition = factor(cell_column,
#                           levels = c("control_stan-119", "within_devi-119", "devi_minus_cont"),
#                           labels = c("control (119ms VOT, observed)", "deviant (119ms VOT, observed)", "deviant /t/ minus control /t/"))) %>%
# filter(condition != "deviant /t/ minus control /t/") %>%
# droplevels()
ylimit <- c(-1.2, 1.2)

fig_diff <- myfunc_permute_diff_plot(df_diff)
fig_diff

### topo ####
fig_topo <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_CHI_identityMMN_English_in_Chinese_blank.jpg")
print(fig_topo)

fig_eng_se <-
  ggdraw() +
  draw_plot(fig_abs, x = 0, y = 2/5, width = 1, height = 3/5) +
  draw_plot(fig_diff, x = 0, y = 0, width = 1, height = 2/5) +
  draw_image(fig_topo, x = 2/3, y = 2/3, width = 1/4, height = 1/5)
fig_eng_se

## English speakers, MC accent ####

sig_times <- df_erp$time[df_erp$time >= 424 & df_erp$time <= 544]
chans <- c('Fz', 'F3', 'F7', 'FC5', 'FC1', 'C3', 'T7', 'CP5', 'CP1', 'Pz', 'P3', 'P7', 'O1', 'Oz', 'P4', 'P8', 'CP6', 'CP2', 'Cz', 'C4', 'FC6', 'FC2', 'F4')

### absolute waves ####
df_abs <- read_delim("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/erp_all.txt", delim = "\t") %>%
  mutate( amp = rowMeans(across(all_of(chans))) ) %>%
  select(participant, group, condition, time, amp) %>%
  filter(condition %in% c("EnCh-devi", "ChEn-stan")) %>%
  mutate(stim_role = case_when(
    condition %in% c("ChEn-stan", "EnCh-stan") ~ "standard",
    condition %in% c("EnCh-devi", "ChEn-devi") ~ "deviant",
    .default = NA)) %>%
  mutate(
    group = factor(group, levels = c("ENG", "CHI", "BILINGUAL"), labels = c("English speakers", "Mandarin-dominant bilinguals", "English-dominant bilinguals")),
    stim_role = factor(stim_role, levels = c("standard", "deviant"))
  ) %>%
  filter(group=="English speakers") %>%
  droplevels()

ylimit <- c(-2, 2)
fig_abs <- myfunc_permute_abs_plot(df_abs, "English bilinguals", "MC accent")
fig_abs

### difference waves ####
df_diff <- df_abs %>%
  select(-stim_role) %>%
  pivot_wider(names_from = condition, values_from = amp) %>%
  mutate(amp = `EnCh-devi` - `ChEn-stan`,
         condition = "deviant - standard")


ylimit <- c(-1.2, 1.2)
fig_diff <- myfunc_permute_diff_plot(df_diff)
fig_diff

### topo ####
fig_topo <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_CHI_identityMMN_English_in_Chinese_blank.jpg")
print(fig_topo)

fig_eng_mc <-
  ggdraw() +
  draw_plot(fig_abs, x = 0, y = 2/5, width = 1, height = 3/5) +
  draw_plot(fig_diff, x = 0, y = 0, width = 1, height = 2/5) +
  draw_image(fig_topo, x = 2/3, y = 2/3, width = 1/4, height = 1/5)

fig_eng_mc

## Mandarin-dominant speakers, SE accent ####

sig_times <- df_erp$time[df_erp$time >= 380 & df_erp$time <= 604]
chans <- c('Fz', 'F3', 'FC5', 'FC1', 'C3', 'CP5', 'CP1', 'Pz', 'P3', 'P7', 'O1', 'O2', 'P4', 'CP6', 'CP2', 'Cz', 'C4', 'FC6', 'FC2', 'F4')

### absolute waves ####
df_abs <- read_delim("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/erp_all.txt", delim = "\t") %>%
  mutate( amp = rowMeans(across(all_of(chans))) ) %>%
  select(participant, group, condition, time, amp) %>%
  filter(condition %in% c("ChEn-devi", "EnCh-stan")) %>%
  mutate(stim_role = case_when(
    condition %in% c("ChEn-stan", "EnCh-stan") ~ "standard",
    condition %in% c("EnCh-devi", "ChEn-devi") ~ "deviant",
    .default = NA)) %>%
  mutate(
    group = factor(group, levels = c("ENG", "CHI", "BILINGUAL"), labels = c("English speakers", "Mandarin-dominant bilinguals", "English-dominant bilinguals")),
    stim_role = factor(stim_role, levels = c("standard", "deviant"))
  ) %>%
  filter(group=="Mandarin-dominant bilinguals") %>%
  droplevels()

ylimit <- c(-2, 2)
fig_abs <- myfunc_permute_abs_plot(df_abs, "Mandarin-dominant bilinguals", "SE accent")
fig_abs

### difference waves ####
df_diff <- df_abs %>%
  select(-stim_role) %>%
  pivot_wider(names_from = condition, values_from = amp) %>%
  mutate(amp = `ChEn-devi` - `EnCh-stan`,
         condition = "deviant - standard")

ylimit <- c(-1.2, 1.2)
fig_diff <- myfunc_permute_diff_plot(df_diff)
fig_diff


fig_topo <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_CHI_identityMMN_English_in_Chinese_blank.jpg")
print(fig_topo)

### topo ####
fig_man_se <-
  ggdraw() +
  draw_plot(fig_abs, x = 0, y = 2/5, width = 1, height = 3/5) +
  draw_plot(fig_diff, x = 0, y = 0, width = 1, height = 2/5) +
  draw_image(fig_topo, x = 2/3, y = 2/3, width = 1/4, height = 1/5)
fig_man_se

# put together figures
fig <-
  ggdraw() +
  draw_plot(fig_eng_se, x = 0, y = 0, width = 1/3, height = 1) +
  draw_plot(fig_eng_mc, x = 1/3, y = 0, width = 1/3, height = 1) +
  draw_plot(fig_man_se, x = 2/3, y = 0, width = 1/3, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), x = c(0, 1/3, 2/3), y = c(1, 1, 1))
fig

ggsave(plot = fig, width = 9, height = 4.5, units = "in", dpi = 300, filename = "OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_wave.jpg")




## permute topos ####
### eng se ####
img_stan <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_ENG_identityMMN_English_in_Chinese_stan.jpg")
print(img_stan)
img_devi <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_ENG_identityMMN_Chinese_in_English_devi.jpg")
print(img_devi)
img_diff <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_ENG_identityMMN_English_in_Chinese_diff.jpg")
print(img_diff)

fig_eng_se <-
  ggdraw() +
  draw_image(img_stan, x = 0, y = 0, width = 1/3, height = 1) +
  draw_image(img_devi, x = 1/3, y = 0, width = 1/3, height = 1) +
  draw_image(img_diff, x = 2/3, y = 0, width = 1/3, height = 1)
fig_eng_se

### eng mc ####
img_stan <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_ENG_identityMMN_Chinese_in_English_stan.jpg")
print(img_stan)
img_devi <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_ENG_identityMMN_Chinese_in_English_devi.jpg")
print(img_devi)
img_diff <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_ENG_identityMMN_Chinese_in_English_diff.jpg")
print(img_stan)

fig_eng_mc <-
  ggdraw() +
  draw_image(img_stan, x = 0, y = 0, width = 1/3, height = 1) +
  draw_image(img_devi, x = 1/3, y = 0, width = 1/3, height = 1) +
  draw_image(img_diff, x = 2/3, y = 0, width = 1/3, height = 1)
fig_eng_mc

### man se ####
img_stan <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_CHI_identityMMN_English_in_Chinese_stan.jpg")
print(img_stan)
img_devi <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_CHI_identityMMN_English_in_Chinese_devi.jpg")
print(img_devi)
img_diff <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo_CHI_identityMMN_English_in_Chinese_diff.jpg")
print(img_stan)

fig_man_se <-
  ggdraw(xlim = c(0, 1), ylim = c(0, 1),clip = "off") +
  draw_image(img_stan, x = 0, y = 0, width = 1/3, height = 1) +
  draw_image(img_devi, x = 1/3, y = 0, width = 1/3, height = 1) +
  draw_image(img_diff, x = 2/3, y = 0, width = 1/3, height = 1)
fig_man_se

# put together
fig <-
  ggdraw() +
  draw_plot(fig_eng_se, x = 0, y = 2/3, width = 1, height = 1/3) +
  draw_plot(fig_eng_mc, x = 0, y = 1/3, width = 1, height = 1/3) +
  draw_plot(fig_man_se, x = 0, y = 0, width = 1, height = 1/3) +
  draw_plot_label(label = c("A", "B", "C"), x = c(0, 0, 0), y = c(1, 2/3, 1/3))
fig

ggsave(plot = fig, width = 6, height = 6, units = "in", dpi = 300, filename = "OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/permute_topo.jpg")

# tfr plots ####
img_ch_stan <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/tfr_ENG_ChEn-stan.jpg")
print(img_ch_stan)
img_en_stan <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/tfr_ENG_EnCh-stan.jpg")
print(img_en_stan)
img_sig_diff <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/tfr_permute_EnglishSpeaker.jpg")
print(img_sig_diff)
img_topo <- image_read("OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/tfr_permute_topo_EnglishSpeaker.jpg")
print(img_topo)

fig <-
  ggdraw() +
  draw_image(img_ch_stan, x = 0, y = 0, width = 1/4, height = 1) +
  draw_image(img_en_stan, x = 1/4, y = 0, width = 1/4, height = 1) +
  draw_image(img_sig_diff, x = 2/4, y = 0, width = 1/4, height = 1) +
  draw_image(img_topo, x = 3/4, y = 0, width = 1/4, height = 1, scale = 0.7) +
  draw_plot_label(label = c("A", "B"), x = c(0, 3/4), y = c(1, 1))
fig

ggsave(plot = fig, width = 6, height = 2, units = "in", dpi = 300, filename = "OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/tfr_permute_oneplot.jpg")

# GAM & brain-language relationship ####

# read in reading data
df_gam <- read.table(file = "~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/gam/df_gam_itemLevel_dependentSE.txt")

df_ppt <- read_delim("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/ppt_demos.txt", delim = "\t") %>%
  filter(exclude != 1) %>%
  select(ppt, CHI, BILINGUAL, ENG, 
         chinese_speaking, chinese_listening, chinese_reading, chinese_writing, chinese_overall, chinese_listeningspeaking,
         english_speaking, english_listening, english_reading, english_writing, english_overall, english_listeningspeaking) %>%
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

df_rating_raw <- read_delim("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/ratings.txt")
df_rating <- na.omit(df_rating_raw) %>%
  inner_join(df_gam_hasPeak, by = "ppt") %>%
  mutate(
    group = factor(group, levels = c("ENG", "CHI", "BILINGUAL"), labels = c("English speakers", "Mandarin-dominant bilinguals", "English-dominant bilinguals")),
    condition = factor(condition, levels = c("ChEn-devi", "EnCh-devi"), labels = c("SE accent", "MC accent"))) %>%
  rename("MC accent familiarity" = familiarity_chinese,
         "MC accent likelihood" = likelihood_chinese,
         "Chinese proficiency" = chinese_overall,
         "English proficiency" = english_overall)

## experience distribution ####
df_dist <- na.omit(df_rating_raw) %>%
  inner_join(df_ppt, by = "ppt") %>%
  mutate(group = case_when(
    ENG==1 ~ "English speakers", # ENG mean English speakers excluding SEA and CBC
    CHI=="1" ~ "Mandarin-dominant bilinguals",
    BILINGUAL=="1" ~ "English-dominant bilinguals",
    .default = "ELSE"
  )) %>%
  filter(group != "ELSE") %>%
  select(ppt, group, familiarity_chinese, likelihood_chinese) %>%
  rename("MC accent familiarity" = familiarity_chinese,
         "MC accent likelihood" = likelihood_chinese) %>%
  pivot_longer(cols = c(
    "MC accent familiarity",
    "MC accent likelihood",
  ), names_to = "exposure_type", values_to = "rating")


df_dist %>% group_by(group, exposure_type) %>%
  summarize(max = max(rating))


## experience correlation plot ####
df_cor_plot <- df_rating %>%
  pivot_longer(cols = c(
    "MC accent familiarity",
    "MC accent likelihood",
  ), names_to = "exposure_type", values_to = "rating") %>%
  # filter(half_area_latency>=150 & half_area_latency<=600) %>%
  # filter(group %in% c("ENG", "BILINGUAL")) %>%
  # filter(NMP < -1) %>%
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
  labs(x = "Rating", y = "Normalized modeled peak") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "top")
print(fig)
# save 
ggsave(plot = fig, width = 7, height = 7, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/cor_NMP_experience.jpg")


## proficiency correlation plot ####
df_cor_plot <- df_rating %>%
  # filter(half_area_latency>=150 & half_area_latency<=600) %>%
  filter(group == "English-dominant bilinguals") %>%
  # filter(NMP < -1) %>%
  droplevels()


fig <-
  ggplot(data = df_cor_plot,
         mapping = aes(x = `Chinese proficiency`, y = trad_erp)) +
  facet_grid(cols = vars(condition)) +
  geom_point(size = 1, alpha = 0.8, color = "#619CFF") +
  # coord_cartesian(ylim = c(-7, 7),
  #                 xlim = c(0, 10)) +
  # scale_x_continuous(breaks = seq(1,10, 1)) +
  geom_smooth(method = "lm",
              se = TRUE,
              linewidth = 0.5,
              alpha = 0.1,
              na.rm = TRUE,
              color = "#619CFF",
              fill = "#619CFF") +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
           label.x.npc = 0, #adjust the label in x axis
           label.y.npc = 0.2, #adjust the label in y axis
           size = 3, color = "#619CFF") +
  labs(x = "English proficiency rating", y = "Traditional ERP (μV)", title = "Mandarin-dominant bilinguals") +
  theme_bw()
print(fig)
# save 
ggsave(plot = fig, width = 7, height = 4.5, units = "in", dpi = 300, filename = "~/OneDrive - University of Toronto/Projects/Yas accent/Manuscript/figures/cor_trad_proficiency.jpg")
