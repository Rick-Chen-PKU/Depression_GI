#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script create graphics diagnostics visualisation for main outcome two-sample MR
# including scatter plot, funnel plot, forest plot and leave-one-out plot

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR")  # your_working_directory

if (!require("pacman")){install.packages("pacman")}
pacman::p_load("TwoSampleMR",  "mr.raps", "MendelianRandomization", "MRPRESSO", "data.table", 
				"stringr", "readxl", "writexl", "stringr", "cowplot")


library(TwoSampleMR)
library(data.table)
library(MendelianRandomization) 
library(writexl)
library(readxl)
library(simex)
library(MRPRESSO)
library(mr.raps)
library(stringr)
library(cowplot)

# sessionInfo('MendelianRandomization')


exposure_trait = 'MDD'
outcome_trait = 'NAFLD'
# 第二部分开始分析
dat = data.frame(fread(paste0('./GWAS/IVS/', exposure_trait, '_to_', outcome_trait, '_instruments.csv')))  

## Scatter plot:
res <- mr(dat, method_list=c('mr_simple_median', 'mr_weighted_median', 'mr_weighted_mode', 'mr_ivw'))
s <- mr_scatter_plot(res, dat)
# #Funnel plot:
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw"))
fu <- mr_funnel_plot(res_single)
# Forest plot
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw"))
fo <- mr_forest_plot(res_single)
## Leave-one-out plot:
res_loo <- mr_leaveoneout(dat)
l <- mr_leaveoneout_plot(res_loo)


# for 2016EA to LCU
new_plot = plot_grid(s[[1]],fu[[1]],fo[[1]],l[[1]],
                     labels = "AUTO", ncol = 2,
                     align = 'h', axis = "l",
                     label_size = 12)
save_plot(paste0("./GWAS/Scatter_plot_",
                 exposure_trait, '_to_', outcome_trait, '.jpg'),
          new_plot,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_height = 6,
          base_width = 6,  # 设置的宽度
          dpi = 800)
