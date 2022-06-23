# neighbourAnalysis - Michael Wozny 20220315
# Process .txt in matlab_neighbour_analysis containing output from MATLAB scripts
#

setwd("~/Documents/Data/neighbour_analysis/alignProject11_ref_16_clean")

# Read in data
first_neighbour_north <- read.delim("./matlab_neighbour_analysis/nearest_neighbour_MOTL_north.txt")
first_neighbour_center <- read.delim("./matlab_neighbour_analysis/nearest_neighbour_MOTL_center.txt")
first_neighbour_south <- read.delim("./matlab_neighbour_analysis/nearest_neighbour_MOTL_south.txt")

# Make directory for concatenated .csv files for strains and summary stats
if (dir.exists("summaryStatistics")==FALSE) dir.create("summaryStatistics")

# Add a second column with the 'neighbourID'
first_neighbour_north[,2] <- 'North'
colnames(first_neighbour_north)[1] <- "Distance"
colnames(first_neighbour_north)[2] <- "NeighbourID"

first_neighbour_center[,2] <- 'Centre'
colnames(first_neighbour_center)[1] <- "Distance"
colnames(first_neighbour_center)[2] <- "NeighbourID"

first_neighbour_south[,2] <- 'South'
colnames(first_neighbour_south)[1] <- "Distance"
colnames(first_neighbour_south)[2] <- "NeighbourID"

my_data <- first_neighbour_north
my_data <- rbind(my_data,first_neighbour_center)
my_data <- rbind(my_data,first_neighbour_south)

# Little prince elephant plots from Michah Allen: https://micahallen.org/2018/03/15/introducing-raincloud-plots/
# Adapted by Michael Wozny 20210411
# For neighbour analysis

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  axis.text.y = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

# Upper and lower SD bounds around the mean
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)
lmad <- function(x) median(x) - mad(x)
umad <- function(x) median(x) + mad(x)

# Summary stats for half-voilin
sum <- ddply(my_data, ~NeighbourID, summarise, mean = mean(Distance), median = median(Distance), mad = mad(Distance), lower = lmad(Distance), upper = umad(Distance), n = length(Distance))
sumFurthest <- ddply(furtherest_neighbour, ~NeighbourID, summarise, mean = mean(Distance), median = median(Distance), mad = mad(Distance), lower = lmad(Distance), upper = umad(Distance), n = length(Distance))

# Order data for x-axis
my_data$NeighbourID <- factor(my_data$NeighbourID, levels = c("South", "Centre", "North"))

g <- ggplot(data = my_data, aes(y = Distance, x = NeighbourID, fill = NeighbourID)) +
  
  # Half-violin of pooled replicates for each strain
  geom_flat_violin(aes(y = Distance, fill = NeighbourID, color = NeighbourID), position = position_nudge(x = .3, y = 0), alpha = .15) +
  
  # Jitter-points of data points coloured by ExperimentDate
  geom_point(aes(y = Distance, color = NeighbourID), position = position_jitterdodge(dodge.width = 1, jitter.height = 0.4), size = 1, alpha = 0.5, shape =16) +
  
  # For black outlines of median points
  geom_point(data = sum, aes(x = NeighbourID, y = median), color = "black", position = position_nudge(x = 0), size = 3.2) +
  
  # Medians and MAD error bars for experimental replicates
  geom_errorbar(data = sum, aes(ymin = lower, ymax = upper, y = median), position = position_nudge(x = 0), width = 0) +
  #geom_point(data = sum, aes(x = NeighbourID, y = median), color = "black", position = position_nudge(x = 0), size = 2.5) +
  geom_point(data = sum, aes(x = NeighbourID, y = median, color = NeighbourID), position = position_nudge(x = 0), size = 2.5) +
  
  # Label number of data-points
  #geom_text(data = sum, aes(x = NeighbourID, label = n, y = 400),position = position_nudge(x = 0),vjust = 0, color = "black") +
  
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values = c("#66c2a5","#8da0cb","#fc8d62")) +
  scale_fill_manual(values = c("#66c2a5","#8da0cb","#fc8d62")) +
  coord_flip() +
  theme_bw() +
  raincloud_theme +
  # Remove x-label 'Strain'
  labs(x = "", y = "Distance to Neighbour (nm)") +
  scale_y_continuous(breaks = get_breaks(by = 10, from = 0),limits = c(0, 80))

# Save as 300 dpi .tff
tiff("neighbours_NCS.tiff", units="in", width=4, height=5, res=300)
g
dev.off()