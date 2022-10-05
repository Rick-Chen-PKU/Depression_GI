#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022

# This script create forestplot for two-sample MR results 


#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("D:/WorkDir/Rstudio/Forestplot")
library(forestplot)
library(dplyr)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("forestplot",  "dplyr")


#---------------------------------------------------------------------#
#                          Forestplot                         	      #----
#---------------------------------------------------------------------#
setwd("D:/WorkDir/Rstudio/Forestplot")
rs_forest_main <- read.csv('forest_main_copy.csv', header = FALSE)
rs_forest_replication <- read.csv('forest_replication_copy.csv', header = FALSE)
rs_forest_combind <- read.csv('forest_combined_copy.csv', header = FALSE)


rs_forest = rs_forest_main
pdf(file = paste0('my_forest_plot', 1,'.pdf'),width = 7,height = 6, onefile = FALSE)
subgps <- setdiff(seq(29),c(1,2,9,16,23))
rs_forest$V1[subgps] <- paste("  ",rs_forest$V1[subgps]) 
forestplot(
  labeltext = as.matrix(rs_forest[,c(1,5,6)]),
  graph.pos = 2, 
  mean = rs_forest$V2, 
  lower = rs_forest$V3,
  upper = rs_forest$V4, 
  is.summary=c(T,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F),
  zero = 1, 
  boxsize = 0.12,
  lineheight = unit(5,'mm'), 
  colgap = unit(3,'mm'),
  lwd.zero = 1.3, 
  lwd.ci = 1, 
  col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
  # col=fpColors(box=c('black'),
  #              summary=c('black'),
  #              lines = c('black'),
  #              zero = c('red')),
  #legend = c("Total effect", "Direct effect", "Indirect effect"),
  #legend_args = fpLegend(pos = list("topright"),title = "Effects",r = unit(.1, "snpc"),gp = gpar(col = "#CCCCCC", lwd = 1.5)),
  xlab="Odds ratio",
  lwd.xaxis=1.2,
  # clip = c(0.50,2),
  xticks = c(0.70, 1.0, 2),
  xlog = TRUE, 
  align=c("l","l"),
  graphwidth = unit(35,"mm"),
  lty.ci = "solid",
  txt_gp=fpTxtGp(label=gpar(cex=0.85), 
                 ticks=gpar(cex=0.7), 
                 xlab=gpar(cex = 0.8)),
  ci.vertices = T
  )
graphics.off()  



rs_forest = rs_forest_replication
pdf(file = paste0('my_forest_plot', 2,'.pdf'),width = 7,height = 6, onefile = FALSE)
subgps <- setdiff(seq(29),c(1,2,9,16,23))
rs_forest$V1[subgps] <- paste("  ",rs_forest$V1[subgps]) 
forestplot(
  labeltext = as.matrix(rs_forest[,c(5,6)]), 
  graph.pos = 1, 
  mean = rs_forest$V2, 
  lower = rs_forest$V3,
  upper = rs_forest$V4, 
  is.summary=c(T,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F),
  zero = 1,
  boxsize = 0.12, 
  lineheight = unit(5,'mm'), 
  colgap = unit(3,'mm'),
  lwd.zero = 1.3,
  lwd.ci = 1, 
  col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
  # col=fpColors(box=c('black'),
  #              summary=c('black'),
  #              lines = c('black'),
  #              zero = c('red')),
  #legend = c("Total effect", "Direct effect", "Indirect effect"),
  #legend_args = fpLegend(pos = list("topright"),title = "Effects",r = unit(.1, "snpc"),gp = gpar(col = "#CCCCCC", lwd = 1.5)),
  xlab="Odds ratio",
  lwd.xaxis=1.2,
  # clip = c(0.50,2),
  xticks = c(0.70, 1.0, 2),
  xlog = TRUE,  
  align=c("l","l"),
  graphwidth = unit(35,"mm"),
  lty.ci = "solid",
  txt_gp=fpTxtGp(label=gpar(cex=0.85), 
                 ticks=gpar(cex=0.7), 
                 xlab=gpar(cex = 0.8)),
  ci.vertices = T
)
graphics.off()  


rs_forest = rs_forest_combind
pdf(file = paste0('my_forest_plot', 3,'.pdf'),width = 7,height = 6, onefile = FALSE)
subgps <- setdiff(seq(29),c(1,2,9,16,23))
rs_forest$V1[subgps] <- paste("  ",rs_forest$V1[subgps]) 
forestplot(
  labeltext = as.matrix(rs_forest[,c(5,6)]), 
  graph.pos = 1, 
  mean = rs_forest$V2, 
  lower = rs_forest$V3,
  upper = rs_forest$V4,
  is.summary=c(T,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F),
  zero = 1, 
  boxsize = 0.12,
  lineheight = unit(5,'mm'),
  colgap = unit(3,'mm'), 
  lwd.zero = 1.3,
  lwd.ci = 1,
  col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
  # col=fpColors(box=c('black'),
  #              summary=c('black'),
  #              lines = c('black'),
  #              zero = c('red')),
  #legend = c("Total effect", "Direct effect", "Indirect effect"),
  #legend_args = fpLegend(pos = list("topright"),title = "Effects",r = unit(.1, "snpc"),gp = gpar(col = "#CCCCCC", lwd = 1.5)),
  xlab="Odds ratio",
  lwd.xaxis=1.2,
  # clip = c(0.50,2),
  xticks = c(0.70, 1.0, 2),
  xlog = TRUE,  
  align=c("l","l"),
  graphwidth = unit(35,"mm"), 
  lty.ci = "solid",
  txt_gp=fpTxtGp(label=gpar(cex=0.85), 
                 ticks=gpar(cex=0.7), 
                 xlab=gpar(cex = 0.8)),
  ci.vertices = T
)
graphics.off()







