# Depression_GI
This code relates to the project titled "Depression and risk of gastrointestinal disorders: a comprehensive two-sample Mendelian randomization study".


# Step 1: Calculate R2 and F-statistic
First, we needed to calculate R2 and F-statistic for MDD instrumental variables. We downloaded complete summary statistics for MDD from PGC.

1_calculate_R2&Fstatistic.R

# Step 2: Find LD proxies for missing SNPs
Next, we needed to find linkage disequilibrium (LD) proxies for genetic variants missing in the outcome dataset. Here, we did this using the LDlinkR package.

2_find_LD_proxies_for_MainConsortium.R

# Step 3: Harmonisation
This step conduct Harmonisation for UK Biobank and FinnGen outcome, produce instrumental variables Output in two-sample MR format.

3.1_harmonisation_for_MainConsortium.R
3.2_harmonisation_for_FinnGen.R

# Step 4: Run two-sample MR analysis
Once our exposure and outcome datasets were ready, we ran the two-sample MR analysis. We did this separately for UK Biobank and FinnGen. 

4_run_2SMR_for_MainConsortium.R
5_run_2SMR_for_FinnGen.R

# Step 5: Run meta-analysis
We combined two-sample MR results obtained in the previous step using a fixed-effect meta-analysis. This was done using the meta package.

5_Mete_analysis.R

# Step 6: Sensitity Analyses
This step focuses on performing sensitivity analysis, including "Remove_pleiotrpy_snps", "CAUSE", "risk_factor_outcome" and "negative_control_outcome".

6_SensitityAnalyses_negative_control_outcome.R
6_SensitityAnalyses_risk_factor_outcome.R
6_SensitityAnalyses_run_2SMR_for_FinnGen_after_remove_pleiotrpy_snps.R
6_SensitityAnalyses_run_2SMR_for_MainConsortium_after_remove_pleiotrpy_snps.R
6_SensitityAnalyses_run_CAUSE_for_MDD_GERD.R
6_SensitityAnalyses_run_CAUSE_for_MDD_IBS.R
6_SensitityAnalyses_run_CAUSE_for_MDD_NAFLD.R
6_SensitityAnalyses_run_CAUSE_for_MDD_PUD.R

# Step 7: Secondary Analyses
This step focuses on performing sensitivity analysis, including "Reverse MR", "LDSC" and "MVMR".

7_SecondaryAnalyses_LDSC.sh
7_SecondaryAnalyses_MVMR.R
7_SecondaryAnalyses_Reversed_direction_find_LD_proxies.R
7_SecondaryAnalyses_Reversed_direction_harmonisation.R
7_SecondaryAnalyses_Reversed_direction_run_2SMR.R

# Step 8: Create visualisations
This step  create graphics diagnostics visualisation for main outcome two-sample MR including scatter plot, funnel plot, forest plot and leave-one-out plot
We created forest plots for discovery, replication and meta analyses using the forestplot packages.

8_Create_graphics_diagnostics_visualisations.R
8_ForestPlot.R

