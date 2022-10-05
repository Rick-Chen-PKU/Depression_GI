#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script conduct Harmonisation for FinnGen, produce IV Output in two-sample MR format

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR/")  # your_working_directory

if (!require("pacman")){install.packages("pacman")}
pacman::p_load("TwoSampleMR",  "mr.raps", "MendelianRandomization", "MRPRESSO", "data.table", 
				"stringr", "readxl", "writexl", "stringr")

#---------------------------------------------------------------------#
#                            Harmonisation                            #----
#---------------------------------------------------------------------#
Project = 'MDD_GI_Project'
path = paste0("/mnt/data1/user/chendongze/Project/", Project, "/disease_")


Harmonisation = function(exposure_list, outcome_trait){ 	
	for (exposure_trait in exposure_list){
		# Reading exposure data
		print(paste0('------Exposure:', exposure_trait, '------'))
		exposure = fread(paste0("./Manual_IVs/", Project, "/", exposure_trait, ".csv")) 
		if ('EAF' %in% colnames(exposure)){
			exp_dat = read_exposure_data(
				filename = paste0("./Manual_IVs/", Project, "/", exposure_trait, ".csv"),
				sep = ",",
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
			exp_dat = read_exposure_data(
				filename = paste0("./Manual_IVs/", Project, "/", exposure_trait, ".csv"),
				sep = ",",
				phenotype_col = "Trait",
				snp_col = "SNP",
				beta_col = "BETA",
				se_col = "SE",
				effect_allele_col = "A1",
				other_allele_col = "A2",
				pval_col = "P",
				samplesize_col = "N")
		}

		
		# Search the IEU OpenGWAS database, fuction "extract_outcome_data" automatically finds the proxies.
		if (outcome_trait == 'GERD'){
			outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'finn-b-K11_REFLUX') 
			outcome_dat$samplesize.outcome = 13141 + 189695  # sample size information can be found in IEU OpenGWAS database
		} else if (outcome_trait == 'IBS'){
			outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'finn-b-K11_IBS')
			outcome_dat$samplesize.outcome = 4605 + 182423  # sample size information can be found in IEU OpenGWAS database
		} else if (outcome_trait == 'PUD'){
			outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'finn-b-K11_GASTRODUOULC')
			outcome_dat$samplesize.outcome = 4510 + 189695  # sample size information can be found in IEU OpenGWAS database
		} else if (outcome_trait == 'NAFLD'){
			outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'finn-b-NAFLD')
			outcome_dat$samplesize.outcome = 894 + 217898  # sample size information can be found in IEU OpenGWAS database
		}
		
		
		dat = TwoSampleMR::harmonise_data(exp_dat, outcome_dat, action=2); nrow(dat)  # defalut action=2 would remove more SNPs, Remove SNPs for being palindromic with intermediate allele frequencies
		# print(paste0('---###palindromic remove SNPs: ', paste(dat[dat$ambiguous == 'TRUE', ]$SNP, collapse=", "), '; total of', length(dat[dat$ambiguous == 'TRUE', ]$SNP)))
		dat = dat[dat$ambiguous == 'FALSE', ]; nrow(dat)  # remove ambiguous SNPs
		dat = TwoSampleMR::steiger_filtering(dat) # steiger_filtering, must contain sample_size column
		# print(paste0('---###steiger_filtering remove SNPs: ', paste(dat[dat$steiger_dir == 'FALSE', ]$SNP, collapse=", "), '; total of', length(dat[dat$steiger_dir == 'FALSE', ]$SNP)))
		dat = dat[dat$steiger_dir == 'TRUE', ]; nrow(dat) # remove SNPs didn't pass steiger filtering test
		write.table(dat, paste0('./Project/',Project, '/', exposure_trait, '_to_', outcome_trait, '_FinnGen_instruments.csv'), col.names=TRUE,row.names = FALSE, sep=",",quote=FALSE)
	}		
}


exposure_list = c('MDD', 'CUD')
outcome_list = c('GERD', 'IBS', 'PUD', 'NAFLD')

for (outcome_trait in outcome_list){
	print(paste0('★★★★★★★Outcome:', outcome_trait, '★★★★★★★'))
	Harmonisation(exposure_list, outcome_trait)
}




