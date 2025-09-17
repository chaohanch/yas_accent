# erp descriptive plot function
myfunc_erp_plot <- function(data, title) {
  ggplot(data, aes(time, amp)) +
  facet_wrap(vars(accent), nrow=2) +
  stat_summary(fun.data = mean_cl_normal, geom = "ribbon", linewidth = 1,
               aes(fill = stim_role),
               alpha = 0.3, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, aes(color = stim_role)) +
  
  # colors and fill
  scale_color_manual(values = c(
    "standard"="#0571b0",
    "deviant"="#ca0020"
  )) +
  scale_fill_manual(values = c(
    "standard"="#0571b0",
    "deviant"="#ca0020"
  )) +
  
  # add lines for x and y zeros
  zero_lines +
  
  # y limits
  coord_cartesian(ylim = ylimit) +
  
  scale_x_continuous(
    breaks = seq(0, 800, by = 200),  # Specify the desired breaks
    # labels = c("(ms)", 0, 200, 400, 600, 800)
  ) +
  scale_y_continuous(
    # labels = number_format(accuracy = 0.1),
    breaks = c(-4, -2, 0, 2, 4),  # Specify the desired breaks
    # labels = c(-1, 0, 1, "(μV)")
  ) +
  labs(
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    # axis.title.x = element_blank(),
    # axis.title.y = element_text(angle = 90, vjust = 0, hjust = -0.1),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.x = element_blank(),
    legend.text = element_blank(),
    # # legend.direction = "vertical",
    legend.position = "none",
    # legend.position.inside = c(0.9, 0.9),
    # legend.margin = margin(0,0,0,0),
    strip.background = element_blank(),
    # strip.text.x = element_text(hjust = 0),
    # text = element_text(size = 16),
    # title = element_text(size = 12),
    # plot.title = element_text(hjust = 0.5,
    #                           margin = margin(b = 0)),
  ) +
  labs(x = "Time (ms)", y = "Amplitude (μV)", title = title)
}


# permutate plot function
myfunc_permute_abs_plot <- function(data, title, subtitle) {
  ggplot(data, aes(time, amp)) +
    stat_summary(fun.data = mean_cl_normal, geom = "ribbon", linewidth = 1,
                 aes(fill = stim_role),
                 alpha = 0.3, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "line", linewidth = 1, aes(color = stim_role)) +
    
    # colors and fill
    scale_color_manual(values = c(
      "standard"="#0571b0",
      "deviant"="#ca0020"
    )) +
    scale_fill_manual(values = c(
      "standard"="#0571b0",
      "deviant"="#ca0020"
    )) +
    
    # add lines for x and y zeros
    zero_lines +
    
    # y limits
    coord_cartesian(ylim = ylimit) +
    
    scale_x_continuous(
      breaks = seq(0, 800, by = 200),  # Specify the desired breaks
      # labels = c("(ms)", 0, 200, 400, 600, 800)
    ) +
    # scale_y_continuous(
    #   # labels = number_format(accuracy = 0.1),
    #   breaks = c(-4, -2, 0, 2, 4),  # Specify the desired breaks
    #   # labels = c(-1, 0, 1, "(μV)")
    # ) +
    labs(y = "Amplitude (μV)", title = title, subtitle = subtitle) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.title.x = element_blank(),
      # axis.title.y = element_text(angle = 90, vjust = 0, hjust = -0.1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_blank(),
      # legend.direction = "vertical",
      legend.position = "none",
      legend.position.inside = c(0.8, 0.8),
      legend.margin = margin(0,0,0,0),
      # text = element_text(size = 16),
      # title = element_text(size = 12),
      # plot.title = element_text(hjust = 0.5,
      #                           margin = margin(b = 0)),
    )
}

myfunc_permute_diff_plot <- function(data) {
  ggplot(data, aes(time, amp)) +
    stat_summary(fun.data = mean_cl_normal, geom = "ribbon", size = 1,
                 fill="black",
                 alpha = 0.3, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "line", linewidth = 1, aes(color = condition)) +
    
    # colors and fill
    scale_color_manual(values = c(
      "black"
    )) +
    
    # add lines for x and y zeros
    zero_lines +
    
    geom_rug(data = data.frame(time = sig_times),
             aes(x = time),
             sides = "b",
             color = "black",
             inherit.aes = FALSE
    ) +
    
    # y limits
    coord_cartesian(ylim = ylimit) +
    
    scale_x_continuous(
      breaks = seq(0, 800, by = 200),  # Specify the desired breaks
      # labels = c("(ms)", 0, 200, 400, 600, 800)
    ) +
    scale_y_continuous(
      # labels = number_format(accuracy = 0.1),
      breaks = c(-1, 0, 1),  # Specify the desired breaks
      # labels = c(-1, 0, 1)
    ) +
    labs(
      x = "Time (ms)", 
      y = "Amplitude (μV)"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      # axis.title.x = element_blank(),
      axis.title.y = element_text(color="white"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      # axis.text.x = element_blank(),
      legend.title = element_blank(),
      # legend.direction = "vertical",
      legend.position = "none",
      legend.position.inside = c(0.8, 0.8),
      legend.margin = margin(0,0,0,0),
      # text = element_text(size = 16),
      # title = element_text(size = 12),
      # plot.title = element_text(hjust = 0.5,
      #                           margin = margin(b = 0)),
    )
}

# add 0 lines ####
zero_lines <- list(
  geom_vline(xintercept = 0, linetype = "solid", alpha = 0.2),
  geom_hline(yintercept = 0, linetype = "solid", alpha = 0.2)
)

