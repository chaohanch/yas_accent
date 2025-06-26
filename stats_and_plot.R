library(tidyverse)
library(Hmisc)
library(ggpubr)
library(rhdf5)
library(lmerTest)

source("flattenCorrMatrix.R")
source("~/Documents/R functions/ggplot_theme_Publication-2.R")
source("~/Documents/R functions/corstars.R")


#### get ERP data ####

# extract channel labels
chan_labels <- read.table("input/channel_names.txt")$V1

# create time labels (1000Hz sampling rate)
time_labels <- read.table("input/time_names.txt")$V1
subj_labels <- read.table("input/subj_names.txt")$V1

# condition labels
condition_labels <- read.table("input/condition_list.txt")$V1

# bad subjects
bad_subjs <- read.table("bad_subjects.txt")$V1

# subject group
english_subjs <- read.table("input/English_subj_list.txt")$V1
chinese_subjs <- read.table("input/Chinese_subj_list.txt")$V1
southasia_subjs <- read.table("input/SouthAsia_subj_list.txt")$V1

#### get data ####

# Open the HDF5 file
erp_data <- h5read("ERP data/erp_data.h5", "erp_data")
# shape the dataset back to (channel x time x cell x subject)
erp_data <- aperm(erp_data, c(4, 3, 2, 1)) * 1e6 # convert volt to microvolt


# get amplitude

#### % Chinese speakers ####

#### %% early cluster ####
# pick electrodes
elecs <- read.table("input/stats_channel_index_CH speakers_ChEn-devi_minus_EnCh-stan_cluster0.txt")$V1
# pick times
times <- read.table("input/stats_time_index_CH speakers_ChEn-devi_minus_EnCh-stan_cluster0.txt")$V1

# pick subjects
pick_subjs <- read.table("input/Chinese_subj_list.txt")$V1
# get subject indices for subjects of interest
subjs <- subj_labels %in% pick_subjs

# average across selected channels and times
tmp <- apply(erp_data[elecs, times, , subjs, drop = FALSE], c(3, 4), mean)

# add labels for each dimension
dimnames(tmp) <- list(
  Dim1 = condition_labels,
  Dim2 = subj_labels[subjs]
)

# convert the data for correlation
df_ch_erp_early <- as.data.frame(t(tmp)) %>%
  # convert row labels to a column
  rownames_to_column(var = "subject") %>%
  # add group info
  mutate(speaker = "Mandarin") %>%
  # compute iMMN
  mutate(
    iMMN_ChEn_early = `ChEn-devi` - `EnCh-stan`,
    iMMN_EnCh_early = `EnCh-devi` - `ChEn-stan`
  ) %>%
  dplyr::select(subject, iMMN_ChEn_early:iMMN_EnCh_early, speaker)

#### %% late cluster ####
# pick electrodes
elecs <- read.table("input/stats_channel_index_CH speakers_ChEn-devi_minus_EnCh-stan_cluster1.txt")$V1
# pick times
times <- read.table("input/stats_time_index_CH speakers_ChEn-devi_minus_EnCh-stan_cluster1.txt")$V1

# pick subjects
pick_subjs <- read.table("input/Chinese_subj_list.txt")$V1
# get subject indices for subjects of interest
subjs <- subj_labels %in% pick_subjs

# average across selected channels and times
tmp <- apply(erp_data[elecs, times, , subjs, drop = FALSE], c(3, 4), mean)

# add labels for each dimension
dimnames(tmp) <- list(
  Dim1 = condition_labels,
  Dim2 = subj_labels[subjs]
)

# convert the data for correlation
df_ch_erp_late <- as.data.frame(t(tmp)) %>%
  # convert row labels to a column
  rownames_to_column(var = "subject") %>%
  # add group info
  mutate(speaker = "Mandarin") %>%
  # compute iMMN
  mutate(
    iMMN_ChEn_late = `ChEn-devi` - `EnCh-stan`,
    iMMN_EnCh_late = `EnCh-devi` - `ChEn-stan`
  ) %>%
  dplyr::select(subject, iMMN_ChEn_late:iMMN_EnCh_late, speaker)

# combine both windows' data
df_ch_erp <- df_ch_erp_early %>%
  inner_join(df_ch_erp_late, by = c("subject", "speaker")) %>%
  mutate(
    iMMN_ChEn = (iMMN_ChEn_early + iMMN_ChEn_late)/2,
    iMMN_EnCh = (iMMN_EnCh_early + iMMN_EnCh_late)/2,
  ) %>%
  dplyr::select(subject, speaker, iMMN_ChEn, iMMN_EnCh)



#### % English speakers ####

# pick electrodes
elecs <- read.table("input/stats_channel_index_EN speakers_ChEn-devi_minus_EnCh-stan_cluster0.txt")$V1
# pick times
times <- read.table("input/stats_time_index_EN speakers_ChEn-devi_minus_EnCh-stan_cluster0.txt")$V1

# pick subjects
pick_subjs <- read.table("input/English_subj_list.txt")$V1
# get subject indices for subjects of interest
subjs <- subj_labels %in% pick_subjs

# average across selected channels and times
tmp <- apply(erp_data[elecs, times, , subjs, drop = FALSE], c(3, 4), mean)

# add labels for each dimension
dimnames(tmp) <- list(
  Dim1 = condition_labels,
  Dim2 = subj_labels[subjs]
)

df_en_erp <- as.data.frame(t(tmp)) %>%
  # convert row labels to a column
  rownames_to_column(var = "subject") %>%
  # add group info
  mutate(speaker = "English") %>%
  # compute iMMN
  mutate(
    iMMN_ChEn = `ChEn-devi` - `EnCh-stan`,
    iMMN_EnCh = `EnCh-devi` - `ChEn-stan`
  ) %>%
  dplyr::select(subject, iMMN_ChEn:iMMN_EnCh, speaker)

# combine the two groups of data
df_erp <- rbind(df_en_erp, df_ch_erp) %>%
  mutate(speaker = as.factor(speaker))


#### % CH+EN speaker ####
# pick electrodes
elecs <- read.table("input/stats_channel_index_EN+CH speakers_All_blocks_cluster0.txt")$V1
# pick times
times <- read.table("input/stats_time_index_EN+CH speakers_All_blocks_cluster0.txt")$V1

# pick subjects
pick_subjs <- c(english_subjs, chinese_subjs)
# get subject indices for subjects of interest
subjs <- (subj_labels %in% pick_subjs) & (!subj_labels %in% bad_subjs)

# average across selected channels and times
tmp <- apply(erp_data[elecs, times, , subjs, drop = FALSE], c(3, 4), mean)

# add labels for each dimension
dimnames(tmp) <- list(
  Dim1 = condition_labels,
  Dim2 = subj_labels[subjs]
)

df_erp <- as.data.frame(t(tmp)) %>%
  # convert row labels to a column
  rownames_to_column(var = "subject") %>%
  # add group info
  mutate(speaker = case_when(subject %in% english_subjs ~ "English",
                             subject %in% chinese_subjs ~ "Mandarin",
                             subject %in% southasia_subjs ~ "SouthAsia",
                             TRUE ~ NA)) %>%
  # compute iMMN
  mutate(
    iMMN_ChEn = `ChEn-devi` - `EnCh-stan`,
    iMMN_EnCh = `EnCh-devi` - `ChEn-stan`,
    iMMN_InEn = `InEn-devi` - `EnIn-stan`,
    iMMN_EnIn = `EnIn-devi` - `EnIn-stan`
  ) %>%
  dplyr::select(subject, iMMN_ChEn:iMMN_EnCh, iMMN_InEn:iMMN_EnIn, speaker)


#### combine data ####

# read in language exposure data
df_language <- read_csv("input/language_exposure_data.csv") %>%
  mutate(subject = as.character(subject))

# combine all data
df_cor <- df_erp %>%
  inner_join(df_language, by = "subject") %>%
  mutate(
  #   en_aoa = case_when(L1 == "English" ~ as.numeric(L1_aoa),
  #                      L2 == "English" ~ as.numeric(L2_aoa),
  #                      L3 == "English" ~ as.numeric(L3_aoa),
  #                      TRUE ~ NA),
    en_use = case_when(L1 == "English" ~ as.numeric(L1_use),
                       L2 == "English" ~ as.numeric(L2_use),
                       L3 == "English" ~ as.numeric(L3_use),
                       TRUE ~ NA),
    ch_use = case_when(L1 == "Mandarin" ~ as.numeric(L1_use),
                       L2 == "Mandarin" ~ as.numeric(L2_use),
                       L3 == "Mandarin" ~ as.numeric(L3_use),
                       TRUE ~ NA),
  #   ch_expo = (Familiarity_Chinese_Accent + Likelihood_Chinese_Accent) / 2,
  ) %>%
  dplyr::select(c(subject:Likelihood_Chinese_Accent, en_use))

# write the csv file
write_csv(file = "df_cor.csv", x = df_cor)

#### correlations ####
# cor
corstars(df_cor[, c(
  "iMMN_ChEn",
  "iMMN_EnCh",
  "Familiarity_Chinese_Accent",
  "Likelihood_Chinese_Accent",
  "en_use"
  # "ch_use"
)])

df_cor_long <- df_cor %>%
  pivot_longer(cols = iMMN_ChEn:iMMN_EnCh, names_to = "stim", values_to = "amp") %>%
  mutate(speaker = factor(speaker, levels = c("Mandarin", "English"), labels = c("Chinese group", "English group")),
         stim = factor(stim, levels = c("iMMN_ChEn", "iMMN_EnCh"), labels = c("SE", "CE")))


#### stats ####
en_ch <- c(-1/2, 1/2)
contrasts(df_cor_long$stim) <- cbind(en_ch)

EN_CH <- c(-1/2, 1/2)
contrasts(df_cor_long$speaker) <- cbind(EN_CH)


model <- lmer(amp ~ speaker * stim + (1|subject), data = df_cor_long)

ggplot(data = df_cor_long) +
  facet_grid(vars(speaker)) +
  geom_histogram(aes(x = amp, fill = stim, color = stim), alpha = 0.5, position = "identity")

#### correlation plotting ####


fig <-
ggplot(data = df_cor_long, mapping = aes(x = stim, y = amp, color = speaker, group = speaker)) +
  geom_point(stat = "summary", fun = mean,
             position = position_dodge(.2)) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2,
                position = position_dodge(.2)) +
  geom_line(stat = "summary", fun = mean,
            position = position_dodge(.2)) +
  scale_colour_Publication() +
  theme_Publication() +
  theme(
    # change font size of non-title/non-label texts
    text = element_text(size = 14),
    # change font size of title/label texts
    title = element_text(size = 12),
    axis.ticks.x = element_blank(), # remove tick marks
    axis.title.x = element_blank(),
    # change the position of legend. "none" to hide legend
    legend.position = c(0.2, 0.9),
    legend.title = element_blank()) +
  # # align legend labels vertically
  guides(color = guide_legend(direction = "vertical")) +
  # coord_cartesian(ylim = c(-3, 3)) +
  labs(y = "Amplitude (μV)")

ggsave(plot = fig, width = 5, height = 4, units = "in", dpi = 300, filename = "Figures/group_effect.png")


#### correlation ####
df_cor_plot <- df_cor_long %>%
  pivot_longer(cols = Familiarity_Chinese_Accent:Likelihood_Chinese_Accent, names_to = "exposure_type", values_to = "rating") %>%
  mutate(exposure_type = factor(exposure_type, levels = c("Familiarity_Chinese_Accent", "Likelihood_Chinese_Accent"), labels = c("Familiarity with\nChinese Accented English", "Likelihood of Exposure to\nChinese Accented English")),
                                rating = as.factor(rating))


fig <-
ggplot(data = na.omit(df_cor_plot),
       mapping = aes(x = rating, y = amp, group = speaker, color = speaker, fill = speaker)) +
  facet_grid(stim ~ exposure_type) +
  geom_point(size = 1, alpha = 0.8) +
  scale_colour_Publication() +
  scale_fill_Publication() +
  theme_Publication() +
  # coord_cartesian(ylim = c(-7, 7),
  #                 xlim = c(0, 10)) +
  # scale_x_continuous(breaks = seq(1,10, 1)) +
  geom_smooth(method = "lm",
              se = TRUE,
              linewidth = 0.5,
              alpha = 0.1,
              na.rm = TRUE) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
           # label.x.npc = 0, #adjust the label in x axis
           # label.y.npc = 1, #adjust the label in y axis
           size = 3) +
  theme(legend.title = element_blank()) +
  labs(x = "Rating (scale: 1-10)", # x-axis label
       y = "Mean iMMN amplitude (μV)") # y-axis label
print(fig)
# save 
ggsave(plot = fig, width = 7, height = 7, units = "in", dpi = 300, filename = "Figures/correlations.png")
