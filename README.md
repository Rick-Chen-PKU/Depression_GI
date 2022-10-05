# Depression_GI
This code relates to the project titled "Depression and risk of gastrointestinal disorders: a comprehensive two-sample Mendelian randomization study".


# Step 1: Create exposure datasets
First, we needed to create exposure datasets. We downloaded complete summary statistics for epigenetic age acceleration (as measured by GrimAge, PhenoAge, HannumAge and Intrinsic HorvathAge) from https://datashare.ed.ac.uk/handle/10283/3645

1_create_exposure_datasets.R

# Step 2: Find LD proxies for missing SNPs
Next, we needed to find linkage disequilibrium (LD) proxies for genetic variants missing in the outcome dataset. Here, we did this using the LDlinkR package.

2.1_find_LD_proxies_for_FinnGen_outcomes.R

2.2_find_LD_proxies_for_PRACTICAL_subtypes.R

2.3_find_LD_proxies_for_CIMBA_subtypes.R


# Step 3: Run two-sample MR analysis
Once our exposure and outcome datasets were ready, we ran the two-sample MR analysis. We did this separately for FinnGen, UK Biobank and international consortiums. We also ran the analysis for cancer subtypes in international consortiums.

3.1_run_2SMR_for_FinnGen.R

3.2_run_2SMR_for_UKB.R

3.3_run_2SMR_for_consortiums.R

3.4_run_2SMR_for_additional_consortium_subtypes.R

# Step 4: Run meta-analysis
We combined two-sample MR results obtained in the previous step using a fixed-effect meta-analysis. This was done using the meta package.

4_run_metaanalysis.R

# Step 5: Create visualisations
We created plots for main, secondary and sensitivity analyses using the meta and ggforestplot packages.

5.1_create_metaanalysis_plots_detail.R

5.2_create_metaanalysis_plots_colour.R

5.3_create_plots_for_subtypes.R

5.4_run_2SMR_for_cancer_byproxy_outcomes_and_create_plots.R

5.5_create_ldsc_plots.R

# Step 6: Run CAUSE analysis
We ran CAUSE analyses for GrimAge acceleration and prostate and colorectal cancers only, as a sensitivity analysis.

6.1_run_CAUSE_for_GrimAge_PrCinPRACTICAL.R

6.2_run_CAUSE_for_GrimAge_CRCinGECCO.R



