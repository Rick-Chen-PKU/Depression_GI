#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script conduct main two-sample MR for negative control outcome
# including IVW, Simple median, Weighted median,Weighted mode, MR-Egger, DIVW, MR-RAPS, MR-PRESSO methods

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR")  # your_working_directory

if (!require("pacman")){install.packages("pacman")}
pacman::p_load("TwoSampleMR",  "mr.raps", "MendelianRandomization", "MRPRESSO", "data.table", 
				"stringr", "readxl", "writexl", "stringr")

#---------------------------------------------------------------------#
#                         Run two-sample MR                           #----
#---------------------------------------------------------------------#
Project = 'MDD_GI_Project'
path = paste0("/mnt/data1/user/chendongze/Project/", Project, "/disease_")


# ieu-a-1001	2016	Years of schooling	SSGAC	293723
# ukb-b-10831	2018	Pack years of smoking	MRC-IEU	142387
# ukb-b-5779	2018	Alcohol intake frequency.	MRC-IEU	462346
# ieu-a-835	2015	Body mass index	GIANT	322154
# ieu-a-61	2015	Waist circumference	GIANT	232101	
# ukb-b-13702	2018	Time spent doing vigorous physical activity	MRC-IEU	64949
# ieu-b-35	2018	C-Reactive protein level	NA	204402


TSMR_RF <- function(exposure_trait, outcome_trait){
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
	if (outcome_trait == 'Years of schooling'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-19560') 
		outcome_dat$samplesize.outcome = 293723  # sample size information can be found in IEU OpenGWAS database
	} else if (outcome_trait == 'Pack years of smoking'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-533')
		outcome_dat$samplesize.outcome = 142387  # sample size information can be found in IEU OpenGWAS database
	} else if (outcome_trait == 'Alcohol intake frequency'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-533')
		outcome_dat$samplesize.outcome = 462346  # sample size information can be found in IEU OpenGWAS database
	} else if (outcome_trait == 'Body mass index'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-533')
		outcome_dat$samplesize.outcome = 322154  # sample size information can be found in IEU OpenGWAS database
	} else if (outcome_trait == 'Waist circumference'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-533')
		outcome_dat$samplesize.outcome = 232101  # sample size information can be found in IEU OpenGWAS database
	} else if (outcome_trait == 'Time spent doing vigorous physical activity'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-533')
		outcome_dat$samplesize.outcome = 64949  # sample size information can be found in IEU OpenGWAS database
	} else if (outcome_trait == 'C-Reactive protein level'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-533')
		outcome_dat$samplesize.outcome = 204402  # sample size information can be found in IEU OpenGWAS database
	}


	dat = TwoSampleMR::harmonise_data(exp_dat, outcome_dat, action=2); nrow(dat)  # defalut action=2 would remove more SNPs, Remove SNPs for being palindromic with intermediate allele frequencies
	dat = dat[dat$ambiguous == 'FALSE', ]; nrow(dat)  # remove ambiguous SNPs
	dat = TwoSampleMR::steiger_filtering(dat) # steiger_filtering, must contain sample_size column
	dat = dat[dat$steiger_dir == 'TRUE', ]; nrow(dat) 
	MRInputObject <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse =  dat$se.outcome, 
	exposure = exposure_trait, outcome = outcome_trait, snps = dat$SNP, effect_allele = dat$effect_allele.exposure, other_allele = dat$other_allele.exposure, eaf = dat$eaf.exposure)
	# IVW
	IVWObject1 <- MendelianRandomization::mr_ivw(MRInputObject,model= "default",robust = FALSE, penalized = FALSE, correl = FALSE, weights ="simple", psi = 0,distribution = "normal",alpha = 0.05)
	# IVWObject1 
	vec1 = c(IVWObject1@Exposure, IVWObject1@Outcome, IVWObject1@SNPs, 'IVW raw', round(IVWObject1@Estimate,3), round(exp(IVWObject1@Estimate),3), 
			paste0("(", as.character(round(IVWObject1@CILower,3)),",", as.character(round(IVWObject1@CIUpper,3)), ")"),
			paste0("(", as.character(round(exp(IVWObject1@CILower),3)),",", as.character(round(exp(IVWObject1@CIUpper),3)), ")"), round(IVWObject1@StdError,3), IVWObject1@Pvalue,  round(IVWObject1@Heter.Stat[1],2), IVWObject1@Heter.Stat[2],
			round((IVWObject1@Heter.Stat[1] - IVWObject1@SNPs + 1 )/IVWObject1@Heter.Stat[1], 3) )
	# IVWObject2 <- MendelianRandomization::mr_ivw(MRInputObject,model= "default",robust = TRUE, penalized = TRUE, correl = FALSE, weights ="simple", psi = 0,distribution = "normal",alpha = 0.05)
	# IVWObject2 
	# vec2 = 	c(IVWObject2@Exposure, IVWObject2@Outcome, IVWObject2@SNPs, 'IVW robust', round(IVWObject2@Estimate,3), round(exp(IVWObject2@Estimate),3),
			# paste0("(", as.character(round(IVWObject2@CILower,3)),",", as.character(round(IVWObject2@CIUpper,3)), ")"),
			# paste0("(", as.character(round(exp(IVWObject2@CILower),3)),",", as.character(round(exp(IVWObject2@CIUpper),3)), ")"), round(IVWObject2@StdError,3), IVWObject2@Pvalue, NA, NA, NA )

	# median-based
	MedianObject1 <-MendelianRandomization::mr_median(MRInputObject,weighting = "simple",distribution ="normal",alpha = 0.05,iterations = 10000,seed = 314159265) 
	# MedianObject1 
	vec2 = 	c(MedianObject1@Exposure, MedianObject1@Outcome, MedianObject1@SNPs, 'Simple median', round(MedianObject1@Estimate,3), round(exp(MedianObject1@Estimate),3),
			paste0("(", as.character(round(MedianObject1@CILower,3)),",", as.character(round(MedianObject1@CIUpper,3)), ")"),
			paste0("(", as.character(round(exp(MedianObject1@CILower),3)),",", as.character(round(exp(MedianObject1@CIUpper),3)), ")"), round(MedianObject1@StdError,3), MedianObject1@Pvalue, NA, NA, NA )
	MedianObject2 <-MendelianRandomization::mr_median(MRInputObject,weighting = "weighted",distribution ="normal",alpha = 0.05,iterations = 10000,seed = 314159265)  
	# MedianObject2
	vec3 = 	c(MedianObject2@Exposure, MedianObject2@Outcome, MedianObject2@SNPs, 'Weighted median', round(MedianObject2@Estimate,3), round(exp(MedianObject2@Estimate),3),
			paste0("(", as.character(round(MedianObject2@CILower,3)),",", as.character(round(MedianObject2@CIUpper,3)), ")"),
			paste0("(", as.character(round(exp(MedianObject2@CILower),3)),",", as.character(round(exp(MedianObject2@CIUpper),3)), ")"), round(MedianObject2@StdError,3), MedianObject2@Pvalue, NA, NA, NA )
	# MedianObject3 <-MendelianRandomization::mr_median(MRInputObject,weighting = "penalized",distribution ="normal",alpha = 0.05,iterations = 10000,seed = 314159265)  
	# MedianObject3

	# Mode-based estimation
	MBEObject <- mr_mbe(MRInputObject, weighting = "weighted", stderror = "delta", phi = 1, seed = 314159265, iterations = 10000, distribution = "normal", alpha = 0.05)
	# MBEObject
	vec4 = 	c(MBEObject@Exposure, MBEObject@Outcome, MBEObject@SNPs, 'Weighted mode', round(MBEObject@Estimate,3), round(exp(MBEObject@Estimate),3),
			paste0("(", as.character(round(MBEObject@CILower,3)),",", as.character(round(MBEObject@CIUpper,3)), ")"), 
			paste0("(", as.character(round(exp(MBEObject@CILower),3)),",", as.character(round(exp(MBEObject@CIUpper),3)), ")"), round(MBEObject@StdError,3), MBEObject@Pvalue, NA, NA, NA )

	# MR-Egger
	EggerObject1 <-mr_egger(MRInputObject,robust = FALSE, penalized = FALSE,correl =FALSE,distribution = "normal",alpha = 0.05)
	# EggerObject1 
	vec5 = 	c(EggerObject1@Exposure, EggerObject1@Outcome, EggerObject1@SNPs, 'MR-Egger', round(EggerObject1@Estimate,3), round(exp(EggerObject1@Estimate),3),
			paste0("(", as.character(round(EggerObject1@CILower.Est,3)),",", as.character(round(EggerObject1@CIUpper.Est,3)), ")"), 
			paste0("(", as.character(round(exp(EggerObject1@CILower.Est),3)),",", as.character(round(exp(EggerObject1@CIUpper.Est),3)), ")"), round(EggerObject1@StdError.Est,3), EggerObject1@Pvalue.Est, NA, NA, NA )
	# EggerObject2 <-mr_egger(MRInputObject,robust = TRUE, penalized = TRUE,correl =FALSE,distribution = "normal",alpha = 0.05) 
	# EggerObject2 

	# Debiased inverse-variance weighted method
	DIVWObject = mr_divw(MRInputObject, over.dispersion = TRUE, alpha = 0.05, diagnostics = FALSE)
	# DIVWObject
	vec6 = 	c(DIVWObject@Exposure, DIVWObject@Outcome, DIVWObject@SNPs, 'DIVW', round(DIVWObject@Estimate,3), round(exp(DIVWObject@Estimate),3),
			paste0("(", as.character(round(DIVWObject@CILower,3)),",", as.character(round(DIVWObject@CIUpper,3)), ")"),
			paste0("(", as.character(round(exp(DIVWObject@CILower),3)),",", as.character(round(exp(DIVWObject@CIUpper),3)), ")"), round(DIVWObject@StdError,3), DIVWObject@Pvalue, NA, NA, NA )

	# MR-RAPS
	res2 = mr.raps.mle(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, over.dispersion=TRUE, diagnostics=FALSE)
	vec7 = c(exposure_trait, outcome_trait, IVWObject1@SNPs, 'MR-RAPS', round(res2$beta.hat,3), round(exp(res2$beta.hat),3), 
				paste0('(', as.character(round(res2$beta.hat-1.96*res2$beta.se, 2)) , ',', as.character(round(res2$beta.hat+1.96*res2$beta.se, 2)), ')'), 
				paste0('(', as.character(round(exp(res2$beta.hat-1.96*res2$beta.se), 2)) , ',', as.character(round(exp(res2$beta.hat+1.96*res2$beta.se), 2)), ')'), 
				round(res2$beta.se,3), as.numeric(res2$beta.p.value), NA, NA, NA)

	# MR-PRESSO
	if (nrow(dat) >= 4){
		res3 = mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
			OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 2000,  SignifThreshold = 0.05, seed=142857)
		vec81 = c(exposure_trait, outcome_trait, IVWObject1@SNPs,'MR-PRESSO:raw', round(res3$`Main MR results`$`Causal Estimate`[1],3), round(exp(res3$`Main MR results`$`Causal Estimate`[1]),2), 
			paste0('(', as.character(round(res3$`Main MR results`$`Causal Estimate`[1]-1.96*res3$`Main MR results`$`Sd`[1], 2)) , ',', as.character(round(res3$`Main MR results`$`Causal Estimate`[1]+1.96*res3$`Main MR results`$`Sd`[1], 2)), ')'), 
			paste0('(', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[1]-1.96*res3$`Main MR results`$`Sd`[1]), 2)) , ',', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[1]+1.96*res3$`Main MR results`$`Sd`[1]), 2)), ')'), 
			round(res3$`Main MR results`$`Sd`[1], 3), res3$`Main MR results`$`P-value`[1], NA, NA, NA)
		vec82 = c(exposure_trait, outcome_trait, IVWObject1@SNPs - length(res3$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),'MR-PRESSO:Outlier-corrected', 
			round(res3$`Main MR results`$`Causal Estimate`[2],3), round(exp(res3$`Main MR results`$`Causal Estimate`[2]),2), 
			paste0('(', as.character(round(res3$`Main MR results`$`Causal Estimate`[2]-1.96*res3$`Main MR results`$`Sd`[2], 2)) , ',', as.character(round(res3$`Main MR results`$`Causal Estimate`[2]+1.96*res3$`Main MR results`$`Sd`[2], 2)), ')'), 
			paste0('(', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[2]-1.96*res3$`Main MR results`$`Sd`[2]), 2)) , ',', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[2]+1.96*res3$`Main MR results`$`Sd`[2]), 2)), ')'), 
			round(res3$`Main MR results`$`Sd`[2], 3), res3$`Main MR results`$`P-value`[2], NA, NA, NA)
	} else {
		vec81 = c(exposure_trait, outcome_trait, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
		vec82 = c(exposure_trait, outcome_trait, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
	}


	out_df = rbind(vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec81, vec82)
	out_df = data.frame(out_df)
	names(out_df) = c('Exposure', 'Outcome', 'N_snp', 'Method', 'beta', 'OR', 'beta_CI', 'OR_CI', 'SE', 'p_value', 'Q_Stat', 'Q_pvalue', 'I_square')

	# sensitivity analysis
	sen2 = mr_pleiotropy_test(dat)  # Horizontal pleiotropy test
	out_df$egger_intercept = sen2$egger_intercept
	out_df$egger_intercept_p_value = sen2$pval


	if (file.exists(paste0('./', Project, '_MR_Result_RiskFactor.xlsx'))){
		old_dat = read_excel(paste0('./', Project, '_MR_Result_RiskFactor.xlsx'))
		new_dat = rbind(old_dat,"",out_df)
		write_xlsx(new_dat, paste0('./', Project, '_MR_Result_RiskFactor.xlsx')) 
	} else {
		write_xlsx(out_df, paste0('./', Project, '_MR_Result_RiskFactor.xlsx')) 
	} 

}


exposure_trait = 'MDD'
# outcome_trait = 'Years of schooling'
# outcome_trait = 'Pack years of smoking'
# outcome_trait = 'Alcohol intake frequency'
# outcome_trait = 'Body mass index'
# outcome_trait = 'Waist circumference'
# outcome_trait = 'Time spent doing vigorous physical activity'
# outcome_trait = 'C-Reactive protein level'
outcome_list = c('Pack years of smoking', 'Alcohol intake frequency', 'Body mass index', 'Waist circumference', 'Time spent doing vigorous physical activity', 'C-Reactive protein level')
for (outcome_trait in outcome_list){
	TSMR_RF(exposure_trait, outcome_trait)
}






