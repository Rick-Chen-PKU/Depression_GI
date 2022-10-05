#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script conduct ran CAUSE analyses for MDD and IBS, as a sensitivity analysis.
# This needs to be run using the terminal

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR")  # your_working_directory

if (!require("pacman")) install.packages("pacman")
if (!require("cause")) devtools::install_github("jean997/cause")   # the development version; rlang＞=1.0.0 is require
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "readr", "ieugwasr", "data.table", "readr")
library(cause)				


#---------------------------------------------------------------------#
#                   Read exposure and outcome                         #----
#---------------------------------------------------------------------#

X1 <- fread('/mnt/data1/user/chendongze/Project/MDD_GI_Project/disease_MDD.txt')  # exposure
X2 <- fread('/mnt/data1/user/chendongze/Project/MDD_GI_Project/disease_IBS.txt') # outcome

X1 = data.frame(X1)
X2 = data.frame(X2)
X1$BETA = log(X1$OR)
X2$BETA = log(X2$OR)


#---------------------------------------------------------------------#
#                            Merge GWAS data                          #----
#---------------------------------------------------------------------#

X <- gwas_merge(X1, X2, snp_name_cols = c("SNP", "SNP"), 
                       beta_hat_cols = c("BETA", "BETA"), 
                       se_cols = c("SE", "SE"), 
                       A1_cols = c("A1", "A1"), 
                       A2_cols = c("A2", "A2"),
					   pval_cols = c("P", "P"))
head(X)

#---------------------------------------------------------------------#
#                    Calculate nuisance parameters                    #----
#---------------------------------------------------------------------#

# only > 100,000 variants
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)   # about 15 min
head(params$mix_grid)

#---------------------------------------------------------------------#
#                                Clump data                           #----
#---------------------------------------------------------------------#

r2_thresh = 0.01
pval_thresh = 1e-3 #should use larger p value, 1e-03 as used for posteriors, not 5e-08
# method1
X_clump <- X %>% rename(rsid = snp, pval = p1) %>% ieugwasr::ld_clump(dat = ., clump_r2 = r2_thresh, clump_p = pval_thresh, 
					plink_bin = genetics.binaRies::get_plink_binary(), pop = "EUR")

top_vars <- X_clump$rsid
print(length(top_vars)) # 1304

# method2 if snp number is larger
# X2 = X2_copy[,c('SNP','CHR')]
# new_X = merge(X, X2, by.x='snp', by.y='SNP')

# new_X1_2 = new_X[new_X$CHR %in% c(1,2), ]
# new_X3_5 = new_X[new_X$CHR %in% seq(3,5), ]
# new_X6_10 = new_X[new_X$CHR %in% seq(6,10), ]
# new_X11_15 = new_X[new_X$CHR %in% seq(11,15), ]
# new_X16_22 = new_X[new_X$CHR %in% seq(16,22), ]

# X_clump1 <- new_X1_2 %>% rename(rsid = snp, pval = p1) %>% ieugwasr::ld_clump(dat = ., clump_r2 = r2_thresh, clump_p = pval_thresh, plink_bin = genetics.binaRies::get_plink_binary(), pop = "EUR")
# X_clump2 <- new_X3_5 %>% rename(rsid = snp, pval = p1) %>% ieugwasr::ld_clump(dat = ., clump_r2 = r2_thresh, clump_p = pval_thresh, plink_bin = genetics.binaRies::get_plink_binary(), pop = "EUR")
# X_clump3 <- new_X6_10 %>% rename(rsid = snp, pval = p1) %>% ieugwasr::ld_clump(dat = ., clump_r2 = r2_thresh, clump_p = pval_thresh, plink_bin = genetics.binaRies::get_plink_binary(), pop = "EUR")
# X_clump4 <- new_X11_15 %>% rename(rsid = snp, pval = p1) %>% ieugwasr::ld_clump(dat = ., clump_r2 = r2_thresh, clump_p = pval_thresh, plink_bin = genetics.binaRies::get_plink_binary(), pop = "EUR")
# X_clump5 <- new_X16_22 %>% rename(rsid = snp, pval = p1) %>% ieugwasr::ld_clump(dat = ., clump_r2 = r2_thresh, clump_p = pval_thresh, plink_bin = genetics.binaRies::get_plink_binary(), pop = "EUR")

# X_clump = rbind(X_clump1,X_clump2,X_clump3,X_clump4,X_clump5)


#---------------------------------------------------------------------#
#                    MR-CAUSE analysis                                #----
#---------------------------------------------------------------------#

# X is unclumped data and variants clumped data
res <- cause(X=X, variants = top_vars, param_ests = params)   # about 4 min
plot(res$sharing)
plot(res$causal)
res$elpd
summary(res, ci_size=0.95)
plot(res)
plot(res, type="data")

png('./cause_plot/CAUSE_MDD_IBS.png', res=300, height=2000, width=3500)
plot(res)
dev.off()

png('./cause_plot/CAUSE_MDD_IBS_data.png', res=300, height=2000, width=3500)
plot(res, type="data")
dev.off()


