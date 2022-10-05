#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script estimates R2 and F-statistics for each of the exposures

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/MDD_GI_Project/') # your_working_directory

# Load required packages
if (!require('pacman')){install.packages('pacman')}
pacman::p_load("data.table", "writexl", "phenoscanner")


#---------------------------------------------------------------------#
#                            Reading Exposure                         #----
#---------------------------------------------------------------------#
# IVs are extracted directly from the source literature and saved in csv format.
MDD_exp_dat <- fread('./MDD.csv')
CUD_exp_dat <- fread('./CUD.csv')


# combine exposures
exp_dat <- rbind(MDD_exp_dat, CUD_exp_dat)

#---------------------------------------------------------------------#
#                          R2 and F-statistic                         #----
#---------------------------------------------------------------------#

# Calculate R2 and F statistics for each exposure dataset
# method 1
exp_dat$MAF <- ifelse(exp_dat$EAF>=0.5, 1-exp_dat$EAF, exp_dat$EAF)
exp_dat$r2 <- 2 * (abs(exp_dat$BETA)) ** 2 * exp_dat$MAF * (1 - exp_dat$MAF) 
exp_dat$F <- exp_dat$r2 * (exp_dat$N - 2) / (1 - exp_dat$r2)

# method 2
# exp_dat$r2 <- (2 * (exp_dat$BETA^2) * exp_dat$EAF * (1 - exp_dat$EAF)) /
  # (2 * (exp_dat$BETA^2) * exp_dat$EAF * (1 - exp_dat$EAF) + 2 * exp_dat$N * exp_dat$EAF * 
	# (1 - exp_dat$EAF) * exp_dat$SE^2)
# exp_dat$F <- exp_dat$r2 * (exp_dat$N - 2) / (1 - exp_dat$r2)

# method 3
# exp_dat$F <- exp_dat$BETA^2 / exp_dat$SE^2
# exp_dat$r2 <- exp_dat$F_stat/(exp_dat$N-2+exp_dat$F_stat)

# Calculate total R2 for each exposure dataset 
r2_func <- function(id)
{
  x <- exp_dat[which(exp_dat$Trait==id),]
  sum(x$r2, na.rm = T)
}

variance_MDD <- r2_func("MDD"); print(variance_MDD) # 1.84%



# Calculate minimum F-statistic for each exposure dataset 
Fmin_func <- function(id)
{
  x <- exp_dat[which(exp_dat$Trait==id),]
  min(x$F, na.rm = T)
}

Fmin_MDD <- Fmin_func("MDD"); print(Fmin_MDD) # 167


# Calculate maximum F-statistic for each exposure dataset 
Fmax_func <- function(id)
{
  x <- exp_dat[which(exp_dat$Trait==id),]
  max(x$F, na.rm = T)
}

Fmax_MDD <- Fmax_func("MDD"); print(Fmax_MDD) # 436

# Calculate median F-statistic for each exposure dataset 
Fmedian_func <- function(id)
{
  x <- exp_dat[which(exp_dat$Trait==id),]
  median(x$F, na.rm = T)
}

Fmedian_MDD <- Fmedian_func("MDD"); print(Fmedian_MDD) # 186


#save
write_xlsx(exp_dat, "exp_data.xlsx")
