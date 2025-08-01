# # Notes
# 
# -   The analysis goal is to extract GAM measures and correlate GAM-based measures with behavioral data.
# -   The model to fit each individual's data is: bam(uV \~ s(time, condition, bs = 'fs', m = 1) + s(time, stim, bs = 'fs', m = 1), data = df, discrete = TRUE)
# -   condition is the interaction term of poa (dorsal vs. glottal), direction (high_to_low vs. low_to_high) x stim_role (standard vs. deviant)

# Libraries

library(tidyverse)
library(mgcv)
library(itsadug)
source("~/Documents/GitHub/yas_accent/my_functions.R")

# Load data and cleaning
# 
# The dataset consists of the following columns:\
# -   group = participant group (vot vs. f0)\
# -   participant = participant number\
# -   time = time in milliseconds (ranges from -200 to 800, in 4 ms bins)\
# -   poa = place of articulation (dorsal vs. glottal)\
# -   stim_role = stimulus role (standard vs. deviant)\
# -   uV = ERP amplitude in microvoltages\
# -   immn_direction = direction for computing iMMN
# 
# # GAM modeling and GAM-based individual measures
# 
# We extract the following GAM-based individual measures:
# 
# -   **trad_erp**: average amplitude of observed data in specified time window
# -   **model_area**: Modelled Area = geometric area (amplitude \* time) under the GAM curve. This measures the area under the peak (or maximum if there is no peak); only looking for positive area's (or negative areas)
# -   **peak_height**: Height Modelled Peak = height of the peak of the GAM smooth, or the highest point if no peak
# -   **NMP**: Normalized Modelled Peak = a measure of robustness of the peak in units of SDs. I.e., how reliably does this subject show the peak? If value is above 1, then the 95% confidence bands do not overlap, and we can be certain the peak is there. If value is between 0 and 1, then there is a lot of variation between the items.
# -   **half_area_latency**: Modelled Area Median Latency = fractional area latency, i.e. latency at 50% of the area (midpoint)
# -   **model_peak_time**: Modelled Peak Latency = latency of the modelled peak
# 
# Additionally we compute the following measures, for extra information:
#   
# -   **trad_norm_erp**: normalized traditional average
# -   **hasPeak**: TRUE/FALSE; is there a peak in the modeled signal in the search window?
# -   **trad_erp**: average of the GAM smooth in the time window
# 
# For all these measures holds that we look at a specified time window or search window. The traditional measure requires a narrower **time window** (e.g. 150-300 ms post stimulus onset), the GAM measures require a **search window** which can be wider (here we use 0 to 800 ms post stimulus onset). We only look for negative peaks.
# 
# The next section of code runs the GAMs for one single participant, extracts the measures and creates the plots (in a separate pdf document).
# 
# The model to fit each individual's data is:
# 
# bam(uV \~ s(time, condition, bs = 'fs', m = 1) + s(time, stim, bs = 'fs', m = 1), discrete = TRUE)


# define search window
search_min = 150; search_max = 800;
# define classic erp window (this can be from permutation test)
trad_min1 <- 200; trad_max1 <- 350
trad_min2 <- 450; trad_max2 <- 600


# collect files
all_files <- list.files("~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/gam/", pattern = "\\.csv$")

# trim and combine files and code participant
chan = c('Cz', 'Fz', 'FC1', 'FC2', 'CP1', 'CP2', 'C3', 'C4', 'Pz')

# list for indexing condition
`ChEn-devi` <- 1; `EnCh-stan` <- 4
`EnCh-devi` <- 3; `ChEn-stan` <- 2
`InEn-devi` <- 7; `EnIn-stan` <- 6
`EnIn-devi` <- 5; `InEn-stan` <- 8

# initialize the full dataframe
df_gam <- data.frame()
for (file in all_files) {
  
  # get participant
  ppt <- strsplit(file, split = "_")[[1]][1]
  
  df <- my_func_get_single_ppt_data(file)
  
  # modeling
  model <- bam(uV ~ 
                 # poa * stim_role +        # fixed effects
                 s(time, by = condition) +  # smooth for each poa x direction x stim_role
                 s(time, item, bs = "fs", m = 1), # random smooth by item
               data = df, 
               discrete = TRUE)  # for large data for speed
  
  # get values and SE for every individual time point
  min_time <- min(model$model[, "time"])
  max_time <- max(model$model[, "time"])
  nval = length(seq(min_time, max_time))
  
  # summary(model)
  
  # %%%%% extract peak height, peak time, and NMP %%%%%%
  
  for (cond in list(c("ChEn-devi", "EnCh-stan"), c("EnCh-devi", "ChEn-stan"), c("InEn-devi", "EnIn-stan"), c("EnIn-devi", "InEn-stan"))) {
    
    tmp_df <- df %>%
      filter(condition %in% cond) %>%
      droplevels() %>%
      group_by(participant, condition, time) %>%
      dplyr::summarize(mean_uV = mean(uV)) %>%
      ungroup() %>%
      pivot_wider(names_from = condition, values_from = mean_uV) %>%
      mutate(diff = get(cond[1]) - get(cond[2]))
    # get traditional erp measures
    trad_erp1 = mean(tmp_df[tmp_df$time>=trad_min1 & tmp_df$time<=trad_max1, ]$diff)
    trad_erp2 = mean(tmp_df[tmp_df$time>=trad_min2 & tmp_df$time<=trad_max2, ]$diff)
    
    # extract modeled standard and deviant data
    devi <- itsadug::get_modelterm(model, select=get(cond[1]), n.grid = nval, as.data.frame = TRUE)
    stan <- itsadug::get_modelterm(model, select=get(cond[2]), n.grid = nval, as.data.frame = TRUE)
    # get difference for MMN
    dat <- devi %>%
      mutate(condition = "mmn",
             fit = devi$fit - stan$fit,
             se.fit = sqrt(devi$se.fit^2 + stan$se.fit^2))
    
    gam_erp1 = mean(dat[dat$time>=trad_min1 & dat$time<=trad_max1, ]$fit)
    gam_erp2 = mean(dat[dat$time>=trad_min2 & dat$time<=trad_max2, ]$fit)
    
    
    # get gam measures
    results <- my_func_extract_gam_measures(dat, search_min, search_max)
    
    # get values
    results$sdat -> sdat
    results$drv -> drv
    results$hasPeak -> hasPeak
    results$area -> area
    results$peak_height -> peak_height
    results$peak_se -> peak_se
    results$NMP -> NMP
    results$peak_time -> peak_time
    results$half_area_latency -> half_area_latency
    results$firsttime -> firsttime
    results$lasttime -> lasttime
    
    
    # get measures for the current condition
    condition <- cond[1]
    tmp_row <- data.frame(ppt, condition, hasPeak, area, peak_height, peak_se, NMP, peak_time, half_area_latency, trad_erp1, trad_erp2, gam_erp1, gam_erp2)
    # add to the extracted data
    df_gam <- rbind(df_gam, tmp_row)
    
    #%%%%% plot for each condition start %%%%%
    # get data for plotting
    dat$fit.plus.se = dat$fit + dat$se.fit
    dat$fit.minus.se = dat$fit - dat$se.fit
    
    # raw ERP data for plotting
    erp_df <- model$model %>%
      group_by(condition, time) %>%
      dplyr::summarize(mean_uV = mean(uV)) %>%
      pivot_wider(names_from = condition, values_from = mean_uV) %>%
      mutate(
        uV_diff = .[[cond[1]]] - .[[cond[2]]]
      )
    
    # plotting
    fig <- 
      ggplot(dat, aes(x = time, y = fit)) +
      # ribbon
      geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit), fill = "skyblue", alpha = 0.3) +
      # gam modelled mmn
      geom_line(color = "blue", linewidth = 1) +
      # shade the search window
      annotate("rect", xmin = search_min, xmax = search_max, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
      # add horizontal and vertical zero lines
      geom_vline(xintercept = 0, linetype = "solid", alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "solid", alpha = 0.2) +
      # add derivative (scaled for visualization)
      geom_line(data = drv, aes(x = time, y = dYdX*100), color = 'red', linetype = "dashed") +
      # raw difference erp
      geom_line(data = erp_df, aes(time, uV_diff), linetype = "dotted", color="black", linewidth=0.5) +
      # title
      labs(title = paste0(ppt, "_", cond[1], "_minus_", cond[2]),
           subtitle = paste0("NMP = ",round(NMP,digits=2)," (uV: ",round(peak_height,digits=2),", SE: ", round(peak_se,digits=2), ")"),
           x = "Time (ms)", y = "uV") +
      theme_bw()
    
    if (!is.na(peak_time)) {
      # df for shades
      segments_df <- sdat %>%
        filter(time>=firsttime & time<=lasttime) %>%
        mutate(x = time, xend = time,
               y = 0, yend = fit)
      fig <- fig +
        # mark peak time
        geom_vline(xintercept=peak_time, linetype = "dashed", linewidth=1) +
        # mark fractional area latency
        geom_vline(xintercept=half_area_latency, linetype = "dotted", linewidth=1) +
        # add shades to the modeled are
        geom_segment(data = segments_df, aes(x=x, xend=xend, y=y, yend=yend), color = "blue", alpha = 0.1)
    }
    
    # print(fig)
    ggsave(plot = fig, width = 8, height = 5, units = "in", dpi = 300, filename = paste0("~/OneDrive - University of Toronto/Projects/Yas accent/figures/gam/", ppt, "_", cond[1], "_minus_", cond[2], ".png"))

  }
}

write.table(df_gam, file = "~/OneDrive - University of Toronto/Projects/Yas accent/data_analysis/gam/df_gam.txt", sep = "\t")
