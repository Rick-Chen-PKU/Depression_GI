#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script finds LD-proxies for main UKB outcomes, which can then be used in two-sample MR analyses
# Output in two-sample MR format 

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/MDD_GI_Project")  # your_working_directory

if (!require("pacman")){install.packages("pacman")}
pacman::p_load("TwoSampleMR", "tidyverse", "dplyr", "data.table",  "LDlinkR", "readxl", "writexl", "stringr")

#---------------------------------------------------------------------#
#                 Read exposure and outcome datasets                  #----
#---------------------------------------------------------------------#

# Read exposure
exp_dat <- read_excel("exp_data.xlsx", col_names = T)

# Format exposure
if ('EAF' %in% colnames(exp_dat)){
		exp_dat = format_data(exp_dat,
			type = "exposure",
			phenotype_col = "Trait",
			snp_col = "SNP",
			eaf_col = "EAF",
			beta_col = "BETA",
			se_col = "SE",
			effect_allele_col = "A1",
			other_allele_col = "A2",
			pval_col = "P",
			samplesize_col = "N")
} else {
		exp_dat = format_data(exp_dat,
			type = "exposure",
			phenotype_col = "Trait",
			snp_col = "SNP",
			beta_col = "BETA",
			se_col = "SE",
			effect_allele_col = "A1",
			other_allele_col = "A2",
			pval_col = "P",
			samplesize_col = "N")
}

# Function to read outcome
out_func_practical <- function(filepath, var)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- fread(paste0(filepath, 'disease_', var, '.txt'), na.strings="")
  outcome_var$BETA <- log(as.numeric(outcome_var$OR))
  outcome_var$Trait = var
  return(outcome_var)
}

filepath = '/mnt/data1/user/chendongze/Project/MDD_GI_Project/'  # your_summary_data_directory
out_gerd <- out_func_practical(filepath, "GERD")
out_ibs <- out_func_practical(filepath, "IBS")
out_pud <- out_func_practical(filepath, "PUD")
out_nafld <- out_func_practical(filepath, "NAFLD")



#---------------------------------------------------------------------#
#              Identify SNPs that need proxies                        #----
#---------------------------------------------------------------------#

# Function to find list of snps in exposure dataset that are missing from the outcome dataset
find_missing_SNP <- function(out_dat) {
    snps_need_proxy <- setdiff(exp_dat$SNP, intersect(exp_dat$SNP, out_dat$SNP))
	return (snps_need_proxy)
}
	
mis_gerd <- find_missing_SNP(out_gerd) 
mis_ibs <- find_missing_SNP(out_ibs) 
mis_pud <- find_missing_SNP(out_pud) 
mis_nafld <- find_missing_SNP(out_nafld) 
 

length(mis_gerd) # check how many snps are in vector
mis_gerd[1] #see snp missing n#1
mis_gerd[2] #see snp missing n#2
mis_gerd[3] #see snp missing n#3

#---------------------------------------------------------------------#
#                    Find proxies for these snps                      #----
#---------------------------------------------------------------------#

# Function to find LD proxy using LDLINK
# We need to make sure the identified proxy SNPs are available in outcome dataset before continuing 
# If SNPs aren't available, you need to find the next best proxy for the missing SNP

my_token = "bf61a241c044"  # Need to request for LDlink API (https://ldlink.nci.nih.gov/?tab=apiaccess)
find_LD_proxy <- function(snps_need_proxy, out_dat) {
		proxy_list = list()
		for (i in 1:length(snps_need_proxy)){
			print(paste0('---### Total ', length(snps_need_proxy), ' SNPs, rank ', i, ' SNP ### ---'))
			temp <- LDproxy(snps_need_proxy[i], "EUR", "r2", token = my_token, file = F, genome_build = "grch37")
			proxy = as.character(temp[temp$R2 >= 0.8,'RS_Number'])
			temp1 <- intersect(proxy, out_dat$SNP)
			if (length(temp1) != 0){
				proxy_list[snps_need_proxy[i]] =  temp1[1]
			} else {
				proxy_list[snps_need_proxy[i]] =  '<<<Nope>>>'
			}
		}
		return (proxy_list)
		# we could change number of proxies we want to find
}

# These proxies imformation would be used in next section !!!
proxy_snp_gerd  <- find_LD_proxy(mis_gerd, out_gerd)  # ('rs3095337'='rs115507122', 'rs138069266'='rs200343480', 'rs13180906'='rs116755193')
proxy_snp_ibs  <- find_LD_proxy(mis_ibs, out_ibs)  # ('rs3095337'='rs115507122', 'rs138069266'='rs200343480', 'rs13180906'='rs116755193')
proxy_snp_pud  <- find_LD_proxy(mis_pud, out_pud)  # ('rs3095337'='rs115507122', 'rs138069266'='rs200343480', 'rs13180906'='rs116755193')
proxy_snp_nafld  <- find_LD_proxy(mis_nafld, out_nafld)  # ('rs3095337'='rs115507122', 'rs138069266'='rs200343480', 'rs13180906'='rs116755193')

# After checking, we found that one SNP "rs17271123" from CUD → PUD can't find proxy. So we finally remove it when analysis!



