## Date Last Modified: 4/2/25
## Author: Adithya Narayanan


## Task: Read and process data from countsdirectory to generate a visual which depicts concurrent outbreak metrics.


rm(list = ls())


# Import Libraries
library(dplyr)
library(ggplot2)


# Set Directory to countsdirectory - change if using different directory instead of default - see pa_simulations.R
directorypath = "/countsdirectory"
setwd(directorypath)

# This script generates a csv file called summarydf which is used to aggregate data across all csv files in countsdirectory.
# These pieces of code remove the file if it exists to allow for re-running.
# Likewise, the plot produced is stored as plot.png and is removed upon re-running. 
# Fresh plot and summary file are produced.
if (file.exists("summary_df.csv")) {
  file.remove("summary_df.csv")
}

if (file.exists("plot.png")) {
  file.remove("plot.png")
}

# Accumulate all csv files into one large data frame and compute means, medians, quantiles etc. for plots
list_csv_files <- list.files(path = directorypath)
df =  do.call(rbind, lapply(list_csv_files, function(x) read.csv(x, stringsAsFactors = FALSE)))


summary_df = df %>% group_by(days_list) %>%
  summarise(iso_mean = mean(iso_today),
            iso_sd = sd(iso_today),
            quar_mean = mean(quar_today),
            quar_sd = sd(quar_today),
            inf_mean = mean(inf_today),
            inf_sd = sd(inf_today), 
            iso_max = max(iso_today),
            iso_min = min(iso_today), 
            quar_max = max(quar_today),
            quar_min = min(quar_today),
            inf_max = max(inf_today),
            inf_min = min(inf_today), 
            iso_q1 = quantile(iso_today, probs = 0.25), 
            iso_q2 = quantile(iso_today, probs = 0.75), 
            quar_q1 = quantile(quar_today, probs = 0.25),
            quar_q2 = quantile(quar_today, probs = 0.75), 
            inf_q1 = quantile(inf_today, probs = 0.25), 
            inf_q2 = quantile(inf_today, probs = 0.75))

summary_df$iso_c1 = summary_df$iso_mean + (1.645/100)*(summary_df$iso_sd)
summary_df$iso_c2 = summary_df$iso_mean - (1.645/100)*(summary_df$iso_sd)
summary_df$quar_c1 = summary_df$quar_mean + (1.645/100)*(summary_df$quar_sd)
summary_df$quar_c2 = summary_df$quar_mean - (1.645/100)*(summary_df$quar_sd)
summary_df$inf_c1 = summary_df$inf_mean + (1.645/100)*(summary_df$inf_sd)
summary_df$inf_c2 = summary_df$inf_mean - (1.645/100)*(summary_df$inf_sd)


# Store Summary df
write.csv(summary_df, "summary_df.csv")


summary_df <- read.csv("summary_df.csv") # Can be used to simply refer to summary df if it already exists without running prior code.

# Produce Plot
plot = ggplot(summary_df, aes(x = days_list, ymin = iso_q2, ymax = iso_q1)) +
  geom_ribbon(fill = '#E44616', alpha = 0.5) +
  geom_line(aes(y = iso_mean, color = "Isolations", size = 0.5), size = 0.5) +
  
  geom_ribbon(aes(x = days_list, ymin = quar_q2, ymax = quar_q1), fill = '#E4B416', alpha = 0.5)+
  geom_line(aes(y = quar_mean, color = "Quarantines", size = 0.5), size = 0.5) +
  
  geom_ribbon(aes(x = days_list, ymin = inf_q2, ymax = inf_q1), fill = '#76B8FC', alpha = 0.5)+
  geom_line(aes(y = inf_mean, color = "Infections", size = 0.5), size = 0.5) +
  
  scale_x_continuous(breaks=pretty(summary_df$day, n = 5)) +
  scale_y_continuous(breaks=pretty(summary_df$quar_q1, n = 1)) +
  scale_color_manual(name= "Legend", values = c('#76B8FC','#E44616','#E4B416')) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(size = 0.5, linetype = 'dashed',
                                    colour = "gray"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme( axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text( color="black",
                                                                                                     size=16.5), axis.text.y = element_text(color="black",
                                                                                                                                            size=16.5)) +
  xlim(0,160) + ylim(0,400) # Adjust Axes here - the second values within xlim, and ylim denote the axes lengths


# Display plot
print(plot)

# Save as PNG file
ggsave("plot.png", w = 3, h = 3, dpi=2000)



