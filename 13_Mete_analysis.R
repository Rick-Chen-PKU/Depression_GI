#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022

# This script runs a meta-analysis of two-sample MR results 
# across UKB, FinnGen and main international consortiums
# This script also applies an BF correction to main meta-analysed IVW results

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("D:/WorkDir/Rstudio/Forestplot")
library(meta)
library(readxl)
library(tibble)
library(writexl)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("meta", "metafor", "openxlsx", "readxl", "writexl", "tidyverse", "dplyr", "tibble", "ggforestplot", "cowplot")

#---------------------------------------------------------------------#
#                   Load two-sample MR results                        #----
#---------------------------------------------------------------------#
dat1 <- read_excel('./MDD_GI_Project_MR_Result_main.xlsx')
dat1$Study <- 'Main'
dat1 <- subset(dat1, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE'))
dat2 <- read_excel('./MDD_GI_Project_MR_Result_FinnGen.xlsx')
dat2$Study <- 'FinnGen'
dat2 <- subset(dat2, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE'))
data <- rbind(dat1, dat2)
data <- data.frame(data)
data <- na.omit(data)


#---------------------------------------------------------------------#
#                           Run meta-analysis                         #----
#---------------------------------------------------------------------#
# Function to run a fixed effect meta-analysis 
meta_func <- function(method_varname, exp_varname, out_varname){
  temp_data <- data[data$Exposure == exp_varname  & data$Outcome==out_varname & data$Method==method_varname,]
  temp_data$beta = as.numeric(temp_data$beta)
  temp_data$SE = as.numeric(temp_data$SE)
  m <- metagen(beta, SE, data=temp_data,studlab=paste(Study), fixed = TRUE, random = FALSE, prediction=TRUE,sm="SMD") 
  #extract values from meta output
  TE.tibble <- as_tibble(m$TE.fixed)
  se.tibble <- as_tibble(m$seTE.fixed)
  p.tibble <- as_tibble(m$pval.fixed)
  #combine tibbles and change column names
  tibble <- cbind(TE.tibble, se.tibble, p.tibble)
  colnames(tibble) <- c("b", "se", "pval")
  #add columns for exposure, outcome and method
  tibble$exposure <- exp_varname
  tibble$outcome <- out_varname
  tibble$method <- method_varname
  return (data.frame(tibble))
}

method_list = c("IVW raw", "Simple median", "Weighted median", "Weighted mode", "MR-Egger", "MR-RAPS", "MR-PRESSO:raw")
exp_list = c('MDD')
out_list = c('GERD','IBS','PUD','NAFLD')
res <- data.frame()
for (method_varname in method_list){
  for (exp_varname in exp_list){
    for (out_varname in out_list){
      tp <- meta_func(method_varname, exp_varname, out_varname)
      res <- rbind(res, tp)
    } 
  }
}



##################################################################
#               Multiple testing correction
##################################################################

# Function to apply multiple testing correction to main meta-analysed IVW results
multiple_testing_func <- function(results_IVW){
  results_IVW$p.fdr<- p.adjust(results_IVW$pval, method = "fdr")
  results_IVW$p.bon<- p.adjust(results_IVW$pval, method = "bonferroni")
  results_IVW$p.hoch <- p.adjust(results_IVW$pval, method = "hochberg")
  results_IVW$p.sig <- ifelse(results_IVW$pval < .05, "*", "")
  results_IVW$p.fdr.sig <- ifelse(results_IVW$p.fdr < .05, "*", "")
  results_IVW$p.bon.sig <- ifelse(results_IVW$p.bon < .05, "*", "")
  results_IVW$p.hoch.sig <- ifelse(results_IVW$p.hoch < .05, "*", "")
  return(results_IVW)
}

IVW_result = res[res$method=='IVW raw', ]
IVW_pval_adj <- multiple_testing_func(IVW_result) 

# Add FDR p-values
res <- merge(res, IVW_pval_adj[, c("p.bon", "method", "outcome")], by=c("method", "outcome"), all.x = T)

# Combine datasets and save meta-analysis results
write_xlsx(res, './meta_analysis_res_incBF.xlsx')


