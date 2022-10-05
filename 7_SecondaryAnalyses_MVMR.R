#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script conduct MVMR

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR")  # your_working_directory

if (!require("pacman")){install.packages("pacman")}
pacman::p_load("TwoSampleMR",  "MVMR", "MendelianRandomization", "MRPRESSO", "data.table", 
				"stringr", "readxl", "writexl", "LDlinkR")

#---------------------------------------------------------------------#
#                            Run MVMR                                 #----
#---------------------------------------------------------------------#


library(TwoSampleMR)
library(MendelianRandomization) 
library(MVMR)
library(data.table)
library(writexl)
library(readxl)
library(stringr)
library(LDlinkR)

# Preparation
pre_path =  "/mnt/data1/user/chendongze/Project/MDD_GI_Project/disease_"
exposure_trait_1 = 'MDD'
exposure_trait_2= 'Insomnia'  # BMI, T2D, Insomnia
outcome_trait = 'NAFLD'  # GERD,PUD,IBS,NAFLD

# Read exposure and outcome dat
exp_1 = fread(paste0(pre_path, exposure_trait_1, '.txt'))
exp_1$BETA = log(exp_1$OR)
exp_1_raw = exp_1
exp_2 = fread(paste0(pre_path, exposure_trait_2, '.txt'))
exp_2$BETA = log(exp_2$OR)
exp_2_raw = exp_2
outcome = fread(paste0(pre_path, outcome_trait, '.txt'))
outcome$BETA = log(outcome$OR)
outcome_raw = outcome

# ReadI IVs
exp_1_ivs = fread(paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/MDD_GI_Project/', exposure_trait_1, '.csv'))
exp_2_ivs = fread(paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/MDD_GI_Project/', exposure_trait_2, '.csv'))

ivs = Reduce(union,  list(exp_1_ivs$SNP, exp_2_ivs$SNP))  # 求多个向量的并集
ivs = data.frame(SNP=ivs); nrow(ivs)
# ivs_independent = ivs
ivs_independent = clump_data(ivs, clump_kb = 10000, clump_r2 = 0.01, pop = "EUR")
nrow(ivs_independent)

exp_1 = exp_1_raw[exp_1_raw$SNP %in% as.character(ivs_independent$SNP), ]
exp_1 = subset(exp_1, select=c('SNP', 'BETA', 'SE', 'P'))
exp_2 = exp_2_raw[exp_2_raw$SNP %in% as.character(ivs_independent$SNP), ]
exp_2 = subset(exp_2, select=c('SNP', 'BETA', 'SE', 'P'))

# exp_2 = extract_outcome_data(snps = as.character(ivs_independent$SNP), outcomes = 'ieu-b-35') 
# exp_2 = subset(exp_2,select=c('SNP','beta.outcome','se.outcome', 'pval.outcome'))
# names(exp_2) = c('SNP', 'BETA', 'SE', 'P')

outcome = outcome_raw[outcome_raw$SNP %in% as.character(ivs_independent$SNP), ]
outcome = subset(outcome, select=c('SNP', 'BETA', 'SE', 'P'))

dat = merge(exp_1, exp_2, by='SNP')
dat = merge(dat, outcome, by='SNP')
nrow(dat)  # common identify in exposure and outcome SNP
dat = dat[dat$P > 5e-8, ]   # 去除outcome p值大于5e-8的位点
nrow(dat)
Num_ivs = nrow(dat)  # 最后使用的工具变量个数

trait_1_specific_num = length(intersect(dat$SNP, exp_1_ivs$SNP))
trait_2_specific_num = length(intersect(dat$SNP, exp_2_ivs$SNP))


rawdat_mvmr1 = subset(dat, select=c('BETA.x','BETA.y', 'SE.x', 'SE.y', 'BETA', 'SE', 'SNP'))
names(rawdat_mvmr1) = c('exp1_beta','exp2_beta', 'exp1_se', 'exp2_se', 'outcome_beta', 'outcome_se', 'SNP')
rawdat_mvmr1 = na.omit(rawdat_mvmr1)

write_xlsx(rawdat_mvmr1, paste0('./Project/MDD_GI_Project/',exposure_trait_1,'-',exposure_trait_2, '_to_', outcome_trait,'_IVs.xlsx'), col_names = TRUE)  # 保存工具变量

F.data <- format_mvmr(BXGs = rawdat_mvmr1[,c(1,2)],
					BYG = rawdat_mvmr1[,5],
					seBXGs = rawdat_mvmr1[,c(3,4)],
					seBYG = rawdat_mvmr1[,6],
					RSID = rawdat_mvmr1[,7])

# head(F.data)
# mvmrcovmatrix<-matrix(c(1, 0.45, 0.45,1), nrow = 2, ncol = 2)  
# Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,3:4])

sres <- strength_mvmr(r_input = F.data, gencov = 0)
# sres <- strength_mvmr(r_input = F.data, gencov = Xcovmat)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
# pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)

res <- ivw_mvmr(r_input = F.data)
# res <- qhet_mvmr(F.data, mvmrcovmatrix, CI = T, iterations = 1000)  # 这一步会花费大量的时间

out_df = cbind(data.frame(confounders_mediators=c(exposure_trait_1,exposure_trait_2), outcome = outcome_trait,  Raw_Trait_NumIVs=c(length(exp_1_ivs$SNP), length(exp_2_ivs$SNP)), 
				Final_Trait_NumIVs=c(trait_1_specific_num, trait_2_specific_num), NumIVs =Num_ivs), as.data.frame(res)[,c(1,2,4)],
				conditional_F_statistic=c(sres$exposure1, sres$exposure2), Qstat=c(pres$Qstat,pres$Qstat), Qpval=c(pres$Qpval, pres$Qpval))


if (file.exists('./MDD_GI_Project_MVMR_Result.xlsx')){
	old_dat = read_excel('./MDD_GI_Project_MVMR_Result.xlsx')
	new_dat = rbind(old_dat, "", out_df)
	write_xlsx(new_dat,'./MDD_GI_Project_MVMR_Result.xlsx') 
} else {
	write_xlsx(out_df,'./MDD_GI_Project_MVMR_Result.xlsx') 
} 
