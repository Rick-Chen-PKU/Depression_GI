#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script conduct LDSC analysis

#!/bin/bash
#@author: Dongze Chen

trait_list=( "MDD" "GERD" "IBS" "PUD" "NAFLD" )
for trait in "${trait_list[@]}"
do
	./munge_sumstats.py --sumstats /mnt/data1/user/chendongze/Project/MDD_GI_Project/disease_${trait}.txt   --chunksize 500000 --out ./sumstats_data/MDD_GI_Project/disease_${trait} --merge-alleles w_hm3.snplist
done


trait_list=( "GERD" "IBS" "PUD" "NAFLD" ) 


for trait in "${trait_list[@]}"
do 
	./ldsc.py --rg ./sumstats_data/MDD_GI_Project/disease_MDD.sumstats.gz,./sumstats_data/MDD_GI_Project/disease_${trait}.sumstats.gz --ref-ld-chr eur_w_ld_chr/  --w-ld-chr eur_w_ld_chr/ --out ./GeneCorrResults/MDD_GI_Project/MDD_${trait}
done