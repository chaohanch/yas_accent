library(tidyverse)
library(plyr)
library(lmerTest)
library(ez)
library(cowplot)
library(ggpubr)
library(ggrepel)
library(gridExtra)

setwd("~/Library/CloudStorage/OneDrive-UniversityofToronto/Projects/Laryngeal/scripts")
folder <- "Mar_25_2025"

# compute slopes from behavioral data ####

# collect files
file_dir <- "../raw data/behavioral/"
list.files(file_dir)->files
files[grepl("csv",files)]->files.list

## bad subjects ####
bad_subjs <- c(
  '0028', # missing block 3 and 4
  '0030', # missing block 3 and 4
  '0034', # missing block 3
  '0040' # missing block 3 and 4
  
  # '0008', # too many bad trials
  # '0024', # too many bad trials
  # '0037', # too many bad trials
  
  # '0005', # non-native
  # '0010', # non-native
  # '0015', # non-native
  # '0020', # non-native
  # '0023', # non-native
  # '0025', # non-native
  # '0026', # left-handed
)

# trim and combine files and code participant
all=NULL
for (i in 1:length(files.list)){
  
  strsplit(files.list[i], "_")[[1]][2] %in% bad_subjs
  
  read.csv(paste0(file_dir, files.list[i]))->temp
  part = gsub("__","_",paste(substr(files.list[i],14,16),substr(files.list[i],9,12),sep="_"))
  temp$part = part
  if(grepl("F0",files.list[i])){
    temp[,c("practice_sounds","category",	"repetition",	"soundstimuli",	"correctAns","prac_key_resp.keys","exptrial_key_resp.keys", "expName",	"part")]->temp
    all = rbind(all, temp)}
  else{
    temp[,c("practice_sounds","category",	"repetition",	"VOTsounds",	"correctAns","prac_key_resp.keys","exptrial_key_resp.keys", "expName",	"part")]->temp
    names(temp)=c("practice_sounds","category",	"repetition",	"soundstimuli",	"correctAns","prac_key_resp.keys","exptrial_key_resp.keys", "expName",	"part")
    all = rbind(all, temp)}
}




# practice (ha vs. a) check -- out of 8 trials, several people got one wrong and one person got two wrong but their behavioral data look fine. So, we'll keep them. 
all$prac_resp = ifelse(all$prac_key_resp.keys%in%c(0,4),"ha","a")
practice = all[all$practice_sounds!="",]
xtabs(~practice_sounds+prac_resp,practice)
practice$target = ifelse(grepl("_ha_",practice$practice_sounds),"ha","a")
practice$acc = ifelse(practice$target == practice$prac_resp, 1, 0)
ddply(practice,.(part),summarize,accuracy = mean(acc))


# main behavioral data
dat = all[all$soundstimuli!="",]
ddply(dat,.(soundstimuli),mutate,
      f0=as.numeric(strsplit(soundstimuli,"_")[[1]][10]),
      vot=as.numeric(strsplit(soundstimuli,"_")[[1]][12]))->dat
dat$Response = ifelse(dat$exptrial_key_resp.keys%in%c("1","q"),"ka", "ga")
dat$resp = ifelse(dat$exptrial_key_resp.keys%in%c("1","q"),1, 0)
dat$cat = ifelse(dat$category=="low_low","ga","ka")
xtabs(~cat+resp+part,dat)

# graphs
## average k response ####
ddply(dat,.(f0,vot,expName),summarize, mean_k=mean(resp))->sum
fig <-
ggplot(sum,aes(vot,f0,fill=mean_k))+geom_tile()+
  scale_fill_gradient(low = "white",high="dark green")+facet_grid(.~expName)
fig
ggsave(plot = fig, width = 4, height = 3, units = "in", dpi = 330, filename = paste0("../figures/", folder, "/behavioral_mean_k_response.png"))


## by-participant k response ####
ddply(dat,.(f0,vot,expName,part),summarize, mean_k=mean(resp))->sum
ggplot(sum,aes(vot,f0,fill=mean_k))+geom_tile()+
  scale_fill_gradient(low = "white",high="dark green")+
  facet_grid(expName~part)

plot_list <- list()
for (subj_ind in 1:length(unique(dat$part))) {
  
  tmp <- ggplot(data = sum[sum$part==unique(sum$part)[subj_ind], ], aes(x = vot, y = f0, fill = mean_k)) +
    geom_tile() +
    scale_fill_gradient(low = "white",high="dark green") +
    coord_cartesian(xlim = c(0, 50),
                    ylim = c(70, 150)) +
    labs(title = unique(dat$part)[subj_ind])
    # facet_grid(expName ~ .)
  
  plot_list[[subj_ind]] <- tmp
  
}

n_f0 <- length(unique(sum[sum$expName=="F0",]$part))

fig <- plot_grid( plotlist = plot_list, ncol = n_f0, nrow = ceiling(length(unique(dat$part))/n_f0) )
fig
ggsave(plot = fig, width = 50, height = 5, units = "in", dpi = 330, limitsize = FALSE,
       filename = paste0("../figures/", folder, "/behavioral_individual_mean_k_response.png"))

## by-participant percent response ####
dat$target = ifelse(dat$category == "low_low","g","k")

plot_list <- list()
for (subj_ind in 1:length(unique(dat$part))) {
  
  tmp <- ggplot(data = dat[dat$part==unique(dat$part)[subj_ind], ], aes(x = target, y = resp)) +
    stat_summary(fun=mean, geom="bar", position=position_dodge()) +
    stat_summary(fun.data = mean_se, geom = "errorbar",position=position_dodge()) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = unique(dat$part)[subj_ind], y = "/k/ response")

  plot_list[[subj_ind]] <- tmp
  
}

n_f0 <- length(unique(dat[dat$expName=="F0",]$part))
fig <- plot_grid( plotlist = plot_list, ncol = n_f0, nrow = ceiling(length(unique(dat$part))/n_f0) )
fig
ggsave(plot = fig, width = 40, height = 5, units = "in", dpi = 330, limitsize = FALSE,
       filename = paste0("../figures/", folder, "/behavioral_individual_k_response_percent.png"))



## violin plot k response ####
dat_summary <- dat %>%
  group_by(expName, part, cat) %>%
  dplyr::summarize(mean_resp = mean(resp)) %>%
  ungroup()
fig <-
ggplot(data = dat_summary, aes(x = cat, y = mean_resp, color = cat)) +
  facet_wrap(~expName) + 
  geom_violin(trim = TRUE, fill = NA) +
  geom_point(size = 2, position = position_dodge(.9), alpha = 0.4) +
  labs(y = "/ka/ response rate",
       x = "Intended voicing") +
  # theme_Publication() +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank()) # change the position of legend. "none" to hide legend
fig
ggsave(plot = fig, width = 4, height = 3, units = "in", dpi = 330, filename = paste0("../figures/", folder, "/behavioral_violin_k_response.png"))

## normalize data ####
dat_norm <- dat %>%
  group_by(expName) %>%
  mutate(f0_n = as.numeric(scale(f0)),
         vot_n = as.numeric(scale(vot)),
         cue_value = ifelse(expName == "F0", f0_n, vot_n))

## plot functions ####
fig <-
  ggplot(dat_norm, aes(cue_value, resp)) + 
  facet_wrap(~expName) +
  # geom_point() + 
  # stat_smooth(se = FALSE)
  stat_smooth(aes(color = expName, fill = expName), method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  scale_color_manual(values = c(
    "F0"="darkblue",
    "VOT"="gold3"
  )) +
  scale_fill_manual(values = c(
    "F0"="darkblue",
    "VOT"="gold3"
  )) +
  # inidividual
  stat_smooth(aes(cue_value, resp, group = part), geom = "line", alpha = 0.7, color = "grey", method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank()) # change the position of legend. "none" to hide legend
fig


## calculate slope for each individual and plot density ####

## Get the subject ID
subjID <- unique(dat$part)

# initialize an summary dataframe
subj_thresh <- data.frame(matrix(ncol = 4, nrow = length(subjID)))
colnames(subj_thresh) <- c("subj", "thresh", "slope", "group")
subj_thresh$subj <- subjID

for (subj_i in subjID) {
  
  # get experiment group
  cue_group <- strsplit(subj_i, split = "_")[[1]][1]
  
  # Extract the data for each subject
  subj_df <- dat_norm %>%
    filter(part == subj_i) %>%
    droplevels()
  # Get the summary data for a single subject
  subj_summary <- subj_df %>%
    group_by(expName, part, cue_value) %>%
    dplyr::summarize(k_perc = mean(resp)) %>%
    ungroup()
  
  ## Run the generalized linear model to get the coefficients for the sigmoid function
  subj_model <- glm(resp ~ cue_value, data = subj_df, family = "binomial")
  
  # Get the intercept for the sigmoid function: 
  b <- as.numeric(subj_model$coefficients[1])
  # Get the slope for the sigmoid function:
  m <- as.numeric(subj_model$coefficients[2])
  # Fit the sigmoid function to the summarized data
  subj_summary$fitted_k_perc <- exp(b + m * subj_summary$cue_value) / (1 + exp(b + m * subj_summary$cue_value))
  # Get the Threshold (the VOT value when Percentage /t/ is 0.5)
  # When P = 0.5, In(P/(1-P)) = 0 = b + mx, so:
  thresh <- -b/m
  # Get the slope of the curve at the Threshold:
  slope <- (exp(b + m * thresh) * m) / (exp(b + m * thresh) + 1)^2
  # # Get the Uncertainty Region (UR), which is the inverse of the slope:
  # UR <- 1/slope_thresX
  
  ## Create a new summarized dataframe for each subject
  subj_thresh[subj_thresh$subj == subj_i, ] <- c(subj_i, thresh, slope, cue_group)
}

# plot slope distribution
fig_dist <-
  ggplot(data = subj_thresh, aes(group = group, fill = group)) +
  # geom_histogram(aes(x = as.numeric(slope_thresh)), binwidth = 0.1, position = "identity", alpha = 0.5, color = "white") +
  geom_density(aes(x = as.numeric(slope)), alpha = 0.5, color = NA) +
  scale_fill_manual(values = c(
    "F0"="darkblue",
    "VOT"="gold3"
  )) +
  scale_x_continuous(limits = c(-1, 3)) +
  # coord_cartesian(xlim = c(-1, 2.5)) +
  theme_bw() +
  theme(
    # change font size of non-title/non-label texts
    # text = element_text(size = 18),
    # change font size of title/label texts
    # title = element_text(size = 18)
  ) #+ 
  # labs(x = "VOT boundary (ms)", # x-axis label
  #      y = "Number of participants",
  #      title = "Distribution of VOT boundaries") # y-axis label
fig_dist


# chisq.test(xtabs(~cat+resp,dat[dat$part=="0001",])) # Italian, sig, but 50% on low-low --> exclude
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0002",])) # English speaking parents, sig 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0003",])) # Hakka, sig
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0004",])) # Bengali, sig
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0005",])) # Mandarin/Cantonese, not sig. 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0006",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0007",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0008",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0009",])) # not sig
# #chisq.test(xtabs(~cat+resp,dat[dat$part=="0010",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0011",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0012",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0013",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0014",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0015",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0016",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0017",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0018",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0019",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0020",])) # exclude - 50% on highF0
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0021",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0022",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0023",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0024",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0025",])) 
# chisq.test(xtabs(~cat+resp,dat[dat$part=="0026",])) 

## cue by-participant random slopes (cues are z-score transformed) ####
subset(dat, dat$expName == "VOT")->VOT
subset(dat, dat$expName == "F0")->F0
VOT$vot.n = scale(VOT$vot)
F0$f0.n = scale(F0$f0)

glmer(resp~vot.n+(1+vot.n|part),VOT,family="binomial")->fit1
summary(fit1)
glmer(resp~f0.n+(1+f0.n|part),F0,family="binomial")->fit2
summary(fit2)
(VOT.slopes=ranef(fit1)$part[2]+fixef(fit1)[2])
(F0.slopes=ranef(fit2)$part[2]+fixef(fit2)[2])
VOT.slopes$subject = row.names(VOT.slopes)
F0.slopes$subject = row.names(F0.slopes)

write.csv(x = VOT.slopes, file = paste0("../analysis data/", folder, "/slopes_vot.csv"))
write.csv(x = F0.slopes, file = paste0("../analysis data/", folder, "/slopes_f0.csv"))

## compare MMN (immn ~ stim x poa, for each group) ####

# read in vot group data
df_vot_wide <- read_csv(paste0("../analysis data/", folder, "/mmn_VOT_0.188_0.36", ".csv"))
# read in F0 group data
df_f0_wide <- read_csv(paste0("../analysis data/", folder, "/mmn_F0_0.188_0.36", ".csv"))

# reorganize data
df_vot <- df_vot_wide %>%
  pivot_longer(cols = withinBlockMMN_glottal_longVOT:identityMMN_dorsal_shortVOT, names_to = "condition", values_to = "amplitude") %>%
  separate(col = condition, into = c("mmn_type", "poa", "stim")) %>%
  mutate(
    poa = factor(poa, levels = c("dorsal", "glottal")),
    stim = factor(stim, levels = c("longVOT", "shortVOT"))
  )

df_f0 <- df_f0_wide %>%
  pivot_longer(cols = withinBlockMMN_glottal_highF0:identityMMN_dorsal_lowF0, names_to = "condition", values_to = "amplitude") %>%
  separate(col = condition, into = c("mmn_type", "poa", "stim")) %>%
  mutate(
    poa = factor(poa, levels = c("dorsal", "glottal")),
    stim = factor(stim, levels = c("highF0", "lowF0"))
  )


# plot and stats

# for each group
for (mmn_group in c("VOT", "F0")) {
  
  # get df name
  dataName <- paste0("df_", tolower(mmn_group))
  
  # for each mmn type
  for (mmn_type in c("identityMMN", "withinBlockMMN")) {
    
    # get data for plot
    df_plot <- get(dataName)[get(dataName)$mmn_type == mmn_type, ]
    
    # contract a data for label
    df_label <- data.frame(
      poa = rep(levels(df_plot$poa), each = 2, times = 1),
      stim = rep(levels(df_plot$stim), each = 1, times = 4)
    ) %>%
      unite(col = "cell", c("poa", "stim"), sep = "_", remove = FALSE) %>%
      # initialize a column for significant marker
      mutate(sig_marker = NA)
    
    # t-test
    for (poa in levels(df_plot$poa)) {
      for (stim in levels(df_plot$stim)) {
        cellName <- paste(mmn_type, poa, stim, sep = "_")
        t_result <- t.test(get(paste0(dataName, "_wide"))[[cellName]], alternative = "less")
        t_marker <- ifelse(t_result$p.value < 0.001, "***",
                            ifelse(t_result$p.value < 0.01, "**",
                                   ifelse(t_result$p.value < 0.05, "*", "ns"))) # 'ns' for non-significant
        df_label$sig_marker <- ifelse(df_label$cell==paste(poa, stim, sep = "_"), t_marker, df_label$sig_marker)
      }
    }
    
    # RM ANOVA 
    mod <- ezANOVA(
      data = df_plot,
      dv = amplitude,
      wid = subject,
      within = .(poa, stim),
      type = 3,
      detailed = FALSE
    )
    print(mod)
    
    # plot and stats
    fig <-
      ggplot(data = df_plot, aes(x = stim, y = amplitude, color = stim)) +
      facet_wrap(~poa) + 
      geom_violin(trim = TRUE, fill = NA) +
      geom_point(size = 2, position = position_dodge(.9), alpha = 0.4) +
      geom_text_repel(aes(label = subject), hjust = 1.2, size = 2) +
      geom_line(aes(group = subject), color = "grey", alpha = 0.5) +
      labs(y = "Amplitude (μV)",
           title = paste(mmn_group, "group,", mmn_type),
           caption = "Markers indicate significance results from a one-sample t-test against 0 for each condition",
           subtitle = paste0(
             "POA (", levels(df_plot$poa)[1], " vs ", levels(df_plot$poa)[2], "): ", "F(", mod$ANOVA$DFn[1], ", ", mod$ANOVA$DFd[1], ") = ", round(mod$ANOVA$F[1],2), ", p = ", round(mod$ANOVA$p[1],3), "\n",
             "STIM (", levels(df_plot$stim)[1], " vs ", levels(df_plot$stim)[2], "): ", "F(", mod$ANOVA$DFn[2], ", ", mod$ANOVA$DFd[2], ") = ", round(mod$ANOVA$F[2],2), ", p = ", round(mod$ANOVA$p[2],3), "\n",
             "POA x STIM: F(", mod$ANOVA$DFn[3], ", ", mod$ANOVA$DFd[3], ") = ", round(mod$ANOVA$F[3],2), ", p = ", round(mod$ANOVA$p[3],3)
           )) +
      geom_text(data = df_label,
                aes(y = max(df_plot$amplitude) + 0.5, label = sig_marker), 
                size = 5, show.legend = FALSE,
                color = "black") +
      # theme_Publication() +
      theme_bw() +
      theme(
        # text = element_text(size = 12),
        # title = element_text(size = 12),
        axis.title.x = element_blank(),
        # legend.position = "none",
        legend.title = element_blank()) # change the position of legend. "none" to hide legend
    # print(fig)
    
    # save
    ggsave(plot = fig, width = 6, height = 5, units = "in", dpi = 330, filename = paste0("../figures/", folder, "/violin_", mmn_group, "_", mmn_type, ".png"))
    
  }
}

## compare MMN (amp ~ role x stim x poa, for each group) ####

# read in vot group data
df_vot_wide <- read.table(paste0("../analysis data/", folder, "/mmn_VOT_0.188_0.36", ".txt"), header = TRUE, sep = ",")
# read in F0 group data
df_f0_wide <- read.table(paste0("../analysis data/", folder, "/mmn_VOT_0.188_0.36", ".txt"), header = TRUE, sep = ",")

# reorganize data
df_vot <- df_vot_wide %>%
  mutate(glottal_stan_shortVOT = glottal_lowStan_highDevi.stan,
         glottal_devi_longVOT = glottal_lowStan_highDevi.devi,
         glottal_stan_longVOT = glottal_highStan_lowDevi.stan,
         glottal_devi_shortVOT = glottal_highStan_lowDevi.devi,
         dorsal_stan_shortVOT = dorsal_lowStan_highDevi.stan,
         dorsal_devi_longVOT = dorsal_lowStan_highDevi.devi,
         dorsal_stan_longVOT = dorsal_highStan_lowDevi.stan,
         dorsal_devi_shortVOT = dorsal_highStan_lowDevi.devi) %>%
  select(subject, glottal_stan_shortVOT:dorsal_devi_shortVOT) %>%
  pivot_longer(cols = glottal_stan_shortVOT:dorsal_devi_shortVOT, names_to = "condition", values_to = "amplitude") %>%
  separate(col = condition, into = c("poa", "role", "stim")) %>%
  mutate(
    poa = factor(poa, levels = c("dorsal", "glottal")),
    role = factor(role, levels = c("stan", "devi")),
    stim = factor(stim, levels = c("shortVOT", "longVOT"))
  )

df_f0 <- df_f0_wide %>%
  mutate(glottal_stan_lowF0 = glottal_lowStan_highDevi.stan,
         glottal_devi_highF0 = glottal_lowStan_highDevi.devi,
         glottal_stan_highF0 = glottal_highStan_lowDevi.stan,
         glottal_devi_lowF0 = glottal_highStan_lowDevi.devi,
         dorsal_stan_lowF0 = dorsal_lowStan_highDevi.stan,
         dorsal_devi_highF0 = dorsal_lowStan_highDevi.devi,
         dorsal_stan_highF0 = dorsal_highStan_lowDevi.stan,
         dorsal_devi_lowF0 = dorsal_highStan_lowDevi.devi) %>%
  select(subject, glottal_stan_lowF0:dorsal_devi_lowF0) %>%
  pivot_longer(cols = glottal_stan_lowF0:dorsal_devi_lowF0, names_to = "condition", values_to = "amplitude") %>%
  separate(col = condition, into = c("poa", "role", "stim")) %>%
  mutate(
    poa = factor(poa, levels = c("dorsal", "glottal")),
    role = factor(role, levels = c("stan", "devi")),
    stim = factor(stim, levels = c("lowF0", "highF0"))
  )


# plot and stats

# for each group
for (mmn_group in c("VOT", "F0")) {
  
  # get df name
  dataName <- paste0("df_", tolower(mmn_group))
  
  # for each mmn type
    
  # get data for plot
  df_plot <- get(dataName)
  
  # contract a data for label
  df_label <- data.frame(
    poa = rep(levels(df_plot$poa), each = 2, times = 1),
    stim = rep(levels(df_plot$stim), each = 1, times = 4)
  ) %>%
    unite(col = "cell", c("poa", "stim"), sep = "_", remove = FALSE) %>%
    # initialize a column for significant marker
    mutate(sig_marker = NA)
  
  # RM ANOVA 
  mod <- ezANOVA(
    data = df_plot,
    dv = amplitude,
    wid = subject,
    within = .(poa, stim, role),
    type = 3,
    detailed = FALSE
  )
  
  # plot and stats
  fig <-
    ggplot(data = df_plot, aes(x = role, y = amplitude, color = stim)) +
    facet_wrap(~poa) + 
    # geom_violin(trim = TRUE, fill = NA) +
    geom_point(stat = "summary", fun = mean, position = position_dodge(.9)) +
    geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, position=position_dodge(.9)) +
    geom_line(stat = "summary", fun = mean, aes(group = stim), position=position_dodge(.9)) +
    labs(y = "Amplitude (μV)",
         title = paste(mmn_group, "group")
        ) +
    theme_bw() +
    theme(
      # text = element_text(size = 12),
      # title = element_text(size = 12),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()) # change the position of legend. "none" to hide legend
    # save
    ggsave(plot = fig, width = 6, height = 5, units = "in", dpi = 330, filename = paste0("../figures/", folder, "/line_", mmn_group, ".png"))
    
}


## mmn ~ slope correlation ####

# combine slope data with MMN data
left_join(VOT.slopes, df_vot_wide)->VOT.coef
left_join(F0.slopes, df_f0_wide)->F0.coef

# substract glottal effects from dorsal effects to estimate the stop category specific effects
VOT.coef$vot_category_effect <- (VOT.coef$identityMMN_dorsal_longVOT + VOT.coef$identityMMN_dorsal_shortVOT)/2 - (VOT.coef$identityMMN_glottal_longVOT + VOT.coef$identityMMN_glottal_shortVOT)/2
F0.coef$f0_category_effect <- (F0.coef$identityMMN_dorsal_highF0 + F0.coef$identityMMN_dorsal_lowF0)/2 - (F0.coef$identityMMN_glottal_highF0 + F0.coef$identityMMN_glottal_lowF0)/2

# correlation plot cue coefficients vs. all MMN stats
plot(VOT.coef)
plot(F0.coef)

all=NULL
# correlation tests
for(i in 3:11){
  temp=cor.test(VOT.coef$vot.n,VOT.coef[,i])
  all=rbind(all,c(names(VOT.coef)[i], round(as.numeric(temp$p.value),4),round(as.numeric(temp$estimate),4)))
}
for(j in 3:11){
  temp=cor.test(F0.coef$f0.n,F0.coef[,j])
  all=rbind(all,c(names(F0.coef)[j], round(as.numeric(temp$p.value),4),round(as.numeric(temp$estimate),4)))
}
all=as.data.frame(all)
names(all)=c("stats","p-value","r-estimate")
write.csv(x = all, file = paste0("../analysis data/", folder, "/correlation_slope_mmn.csv"))

# converet o long format
vot_coef_long <- VOT.coef %>%
  select(-c(withinBlockMMN_glottal_longVOT:withinBlockMMN_dorsal_shortVOT)) %>%
  pivot_longer(cols = identityMMN_glottal_longVOT:identityMMN_dorsal_shortVOT, names_to = "condition", values_to = "amplitude") %>%
  separate(col = condition, into = c("mmn_type", "poa", "stim")) %>%
  mutate(
    poa = factor(poa, levels = c("dorsal", "glottal")),
    stim = factor(stim, levels = c("longVOT", "shortVOT")),
    group = "vot"
  ) %>%
  dplyr::rename("slope" = vot.n,
                "category_effect" = vot_category_effect)

f0_coef_long <- F0.coef %>%
  select(-c(withinBlockMMN_glottal_highF0:withinBlockMMN_dorsal_lowF0)) %>%
  pivot_longer(cols = identityMMN_glottal_highF0:identityMMN_dorsal_lowF0, names_to = "condition", values_to = "amplitude") %>%
  separate(col = condition, into = c("mmn_type", "poa", "stim")) %>%
  mutate(
    poa = factor(poa, levels = c("dorsal", "glottal")),
    stim = factor(stim, levels = c("highF0", "lowF0")),
    group = "f0"
  ) %>%
  dplyr::rename(slope = "f0.n",
                "category_effect" = f0_category_effect)

df_corPlot <- rbind(vot_coef_long, f0_coef_long) %>%
  mutate(poa = as.factor(poa),
         group = as.factor(group),
         stim = as.factor(stim)) %>%
  na.omit()

mmn_group <- "vot"
fig_corVOT <-
  ggplot(data = df_corPlot[df_corPlot$group==mmn_group,],
         mapping = aes(x = slope, y = amplitude,
                       group = stim,
                       color = stim,
                       fill = stim,
         )) +
  facet_wrap(~poa) +
  geom_point(size = 2, na.rm = TRUE) +
  # geom_text_repel(aes(label = subject), hjust = 1.2, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              alpha = 0.1,
              linewidth = 0.7,
              na.rm = TRUE) +
  # ylim(-0.1, 0.4) +
  # coord_cartesian(ylim = c(-0.1, 0.4)) +
  stat_cor(size = 3, na.rm = TRUE) +
  labs(x = "Slope", # x-axis label
       y = "Amplitude (μV)",
       title = paste(mmn_group, "group")) + # y-axis label
  # geom_text(data = df_txt, aes(x = ewm_x, y = ewm_y, label = ewm_labels, color = Group), size = 4, hjust = 0) +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) # change the position of legend. "none" to hide legend
fig_corVOT

mmn_group <- "f0"
fig_corF0 <-
  ggplot(data = df_corPlot[df_corPlot$group==mmn_group,],
         mapping = aes(x = slope, y = amplitude,
                       group = stim,
                       color = stim,
                       fill = stim,
         )) +
  facet_wrap(~poa) +
  geom_point(size = 2, na.rm = TRUE) +
  # geom_text_repel(aes(label = subject), hjust = 1.2, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              alpha = 0.1,
              linewidth = 0.7,
              na.rm = TRUE) +
  stat_cor(size = 3, na.rm = TRUE) +
  # ylim(-0.1, 0.4) +
  # coord_cartesian(ylim = c(-0.1, 0.4)) +
  labs(x = "Slope", # x-axis label
       y = "Amplitude (μV)",
       title = paste(mmn_group, "group")) + # y-axis label
  # geom_text(data = df_txt, aes(x = ewm_x, y = ewm_y, label = ewm_labels, color = Group), size = 4, hjust = 0) +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) # change the position of legend. "none" to hide legend
fig_corF0

# combine figures
fig <-
  ggdraw() +
  draw_plot(fig_corVOT, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(fig_corF0, x = 0.5, y = 0, width = 0.5, height = 1)
print(fig)

# save
ggsave(plot = fig, width = 10, height = 4, units = "in", dpi = 330, filename = paste0("../figures/", folder, "/correlation_slope_mmn.png"))



# correlation on category effect ####
# converet o long format
vot_coef_long <- VOT.coef %>%
  select(c(vot.n, subject, vot_category_effect)) %>%
  mutate(
    group = "vot"
  ) %>%
  dplyr::rename("slope" = vot.n, "category_effect" = vot_category_effect)

f0_coef_long <- F0.coef %>%
  select(c(f0.n, subject, f0_category_effect)) %>%
  mutate(
    group = "f0"
  ) %>%
  dplyr::rename("slope" = f0.n, "category_effect" = f0_category_effect)

df_corPlot <- rbind(vot_coef_long, f0_coef_long) %>%
  mutate(group = as.factor(group)) %>%
  na.omit()

fig <-
  ggplot(data = df_corPlot,
         mapping = aes(x = slope, y = category_effect,
                       group = group,
                       color = group,
                       fill = group,
         )) +
  geom_point(size = 2, na.rm = TRUE) +
  # geom_text_repel(aes(label = subject), hjust = 1.2, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              alpha = 0.1,
              linewidth = 0.7,
              na.rm = TRUE) +
  # ylim(-0.1, 0.4) +
  # coord_cartesian(ylim = c(-0.1, 0.4)) +
  stat_cor(size = 3, na.rm = TRUE) +
  labs(x = "Slope", # x-axis label
       y = "Amplitude (μV)",
       title = "category effect") + # y-axis label
  # geom_text(data = df_txt, aes(x = ewm_x, y = ewm_y, label = ewm_labels, color = Group), size = 4, hjust = 0) +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) # change the position of legend. "none" to hide legend

print(fig)

# save
ggsave(plot = fig, width = 4, height = 4, units = "in", dpi = 330, filename = paste0("../figures/", folder, "/correlation_category_effect.png"))


# correlation slope ~ within-block mmn, collapsing mmn direction ####
# converet o long format
vot_coef_long <- VOT.coef %>%
  mutate(
    group = "vot",
    dorsal_mmn = (identityMMN_dorsal_longVOT + identityMMN_dorsal_shortVOT)/2,
    glottal_mmn = (identityMMN_glottal_shortVOT + identityMMN_glottal_longVOT)/2,
  ) %>%
  select(c(group, vot.n, subject, dorsal_mmn, glottal_mmn)) %>%
  dplyr::rename("slope" = vot.n)

f0_coef_long <- F0.coef %>%
  mutate(
    group = "f0",
    dorsal_mmn = (identityMMN_dorsal_highF0 + identityMMN_dorsal_lowF0)/2,
    glottal_mmn = (identityMMN_glottal_highF0 + identityMMN_glottal_lowF0)/2,
  ) %>%
  select(c(group, f0.n, subject, dorsal_mmn, glottal_mmn)) %>%
  dplyr::rename("slope" = f0.n)

df_corPlot <- rbind(vot_coef_long, f0_coef_long) %>%
  mutate(group = as.factor(group)) %>%
  na.omit()

fig <-
  ggplot(data = df_corPlot,
         mapping = aes(x = slope, y = dorsal_mmn,
                       group = group,
                       color = group,
                       fill = group,
         )) +
  geom_point(size = 2, na.rm = TRUE) +
  # geom_text_repel(aes(label = subject), hjust = 1.2, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              alpha = 0.1,
              linewidth = 0.7,
              na.rm = TRUE) +
  # ylim(-0.1, 0.4) +
  # coord_cartesian(ylim = c(-0.1, 0.4)) +
  stat_cor(size = 3, na.rm = TRUE) +
  labs(x = "Slope", # x-axis label
       y = "Amplitude (μV)",
       title = "category effect") + # y-axis label
  # geom_text(data = df_txt, aes(x = ewm_x, y = ewm_y, label = ewm_labels, color = Group), size = 4, hjust = 0) +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) # change the position of legend. "none" to hide legend

print(fig)

fig <-
  ggplot(data = df_corPlot,
         mapping = aes(x = slope, y = glottal_mmn,
                       group = group,
                       color = group,
                       fill = group,
         )) +
  geom_point(size = 2, na.rm = TRUE) +
  # geom_text_repel(aes(label = subject), hjust = 1.2, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              alpha = 0.1,
              linewidth = 0.7,
              na.rm = TRUE) +
  # ylim(-0.1, 0.4) +
  # coord_cartesian(ylim = c(-0.1, 0.4)) +
  stat_cor(size = 3, na.rm = TRUE) +
  labs(x = "Slope", # x-axis label
       y = "Amplitude (μV)",
       title = "category effect") + # y-axis label
  # geom_text(data = df_txt, aes(x = ewm_x, y = ewm_y, label = ewm_labels, color = Group), size = 4, hjust = 0) +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) # change the position of legend. "none" to hide legend

print(fig)


# correlation slope ~ immn, collapsing mmn direction ####
# converet o long format
vot_coef_long <- VOT.coef %>%
  mutate(
    group = "vot",
    dorsal_mmn = (withinBlockMMN_dorsal_longVOT + withinBlockMMN_dorsal_shortVOT)/2,
    glottal_mmn = (withinBlockMMN_glottal_shortVOT + withinBlockMMN_glottal_longVOT)/2,
  ) %>%
  select(c(group, vot.n, subject, dorsal_mmn, glottal_mmn)) %>%
  dplyr::rename("slope" = vot.n)

f0_coef_long <- F0.coef %>%
  mutate(
    group = "f0",
    dorsal_mmn = (withinBlockMMN_dorsal_highF0 + withinBlockMMN_dorsal_lowF0)/2,
    glottal_mmn = (withinBlockMMN_glottal_highF0 + withinBlockMMN_glottal_lowF0)/2,
  ) %>%
  select(c(group, f0.n, subject, dorsal_mmn, glottal_mmn)) %>%
  dplyr::rename("slope" = f0.n)

df_corPlot <- rbind(vot_coef_long, f0_coef_long) %>%
  mutate(group = as.factor(group)) %>%
  na.omit()

fig <-
  ggplot(data = df_corPlot,
         mapping = aes(x = slope, y = dorsal_mmn,
                       group = group,
                       color = group,
                       fill = group,
         )) +
  geom_point(size = 2, na.rm = TRUE) +
  # geom_text_repel(aes(label = subject), hjust = 1.2, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              alpha = 0.1,
              linewidth = 0.7,
              na.rm = TRUE) +
  # ylim(-0.1, 0.4) +
  # coord_cartesian(ylim = c(-0.1, 0.4)) +
  stat_cor(size = 3, na.rm = TRUE) +
  labs(x = "Slope", # x-axis label
       y = "Amplitude (μV)",
       title = "dorsal") + # y-axis label
  # geom_text(data = df_txt, aes(x = ewm_x, y = ewm_y, label = ewm_labels, color = Group), size = 4, hjust = 0) +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) # change the position of legend. "none" to hide legend

print(fig)

fig <-
  ggplot(data = df_corPlot,
         mapping = aes(x = slope, y = glottal_mmn,
                       group = group,
                       color = group,
                       fill = group,
         )) +
  geom_point(size = 2, na.rm = TRUE) +
  # geom_text_repel(aes(label = subject), hjust = 1.2, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              alpha = 0.1,
              linewidth = 0.7,
              na.rm = TRUE) +
  # ylim(-0.1, 0.4) +
  # coord_cartesian(ylim = c(-0.1, 0.4)) +
  stat_cor(size = 3, na.rm = TRUE) +
  labs(x = "Slope", # x-axis label
       y = "Amplitude (μV)",
       title = "glottal") + # y-axis label
  # geom_text(data = df_txt, aes(x = ewm_x, y = ewm_y, label = ewm_labels, color = Group), size = 4, hjust = 0) +
  theme_bw() +
  theme(
    # text = element_text(size = 12),
    # title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) # change the position of legend. "none" to hide legend

print(fig)
