#!/bin/sh -e
#######################################################################################################
# Scripts provided by discosnprad, pipeline made by Ben Gibbons.
#######################################################################################################

# Full filtering of a clustered vcf file using various discosnprad scripts and bcftools to arrive at fully filtered vcf file. 

# discosnprad scripts available here: https://github.com/GATB/DiscoSnp/tree/master/discoSnpRAD/post-processing_scripts.
# Clustered vcf files used in this pipeline were already clustered by Rob Elshire using the discosnprad clustering scripts available here: https://github.com/GATB/DiscoSnp/tree/master/discoSnpRAD/clustering_scripts.

# This script is an example of how I filtered my data set for use in down stream analyses for use in distance based and frequency based analyses. 


# First filtered by individual cover -c (minimum depth to be included), Max % missing data -m and Minor allele frequency -f (maf).

python3 filter_vcf_by_indiv_cov_max_missing_and_maf.py -i discoRad_k_31_c_3_D_0_P_5_m_5_clustered.vcf -o ecology_c_8_f_0.05_m_0.vcf -c 8 -f 0.05 -m 0 -s 1
#     OUTPUT
# filter_vcf_by_indiv_cov_max_missing_and_maf.py
# input_file : discoRad_k_31_c_3_D_0_P_5_m_5_clustered.vcf
# output_file : ecology_c_8_f_0.05_m_0.vcf
# filter parameters : indiv_DP>=8, missing<=0 (prop*n_geno = 0.0 * 83), snp-only=1
# 2545929 seen variants, 24785 variants after filtering
# initial missing count (on snps only if snp-only=1) = 148782220, genotypes changed to NA = 25607726, final missing count in new vcf = 0
# initial missing percent (on snps only if snp-only=1) = 70.41 %, final missing percent in new vcf = 0.00 %

# Then filter the file by the single snppercluster script (to account for SNP linkage).
python3 SSPC_script.py -i ecology_c_8_f_0.05_m_0.vcf -o ecology_dp8_sspc.vcf

# Then the heterozygozity paralog filter that exludes heterozygous sites that are at 75% of higher.
python3 filter_paralogs.py.py -i ecology_dp8_sspc.vcf -o ecology_dp8_sspc_paralog_x0.75.vcf -x 0.75 -y 1
#    OUTPUT
# 11964 on 13011 clusters had less than 100.0% of SNP with less than 75.0% heterygous genotypes

# Now to use bcftools to identify 'bad SNP sites' from a triplicate included in the data set, to add a maximum depth cutoff and to remove unwanted samples. 

conda activate bcftools 

# To be able to identify 'bad SNP sites' the vcf file needs to be bgzip'd and then tabix'd to extract these bad sites.
bgzip ecology_dp8_sspc_paralog_x0.75.vcf

tabix -p vcf ecology_dp8_sspc_paralog_x0.75.vcf.gz

# To identify 'bad SNP sites' that may be present. Sample_file.txt conatins the name of the triplicate samples.
bcftools view -S sample_file.txt ecology_dp8_sspc_paralog_x0.75.vcf.gz |  bcftools query -i 'N_PASS(GT="AA")<3 & N_PASS(GT="RR")<3 & N_PASS(GT="het")<3' -f '%CHROM\t%POS\n' > bad_sites2.tsv

# To see how many 'bad SNP sites' there are in the file.
wc -l bad_sites2.tsv
2506 bad_sites2.tsv

#Excluding 'bad SNP sites' from the fully filtered VCF file.
bcftools view -T ^bad_sites2.tsv -O z -o ecology_dp8_sspc_paralog_x0.75_nobadSites.vcf.gz ecology_dp8_sspc_paralog_x0.75.vcf.gz

# Formating filtered vcf file to have a lower and upper limit for depth (this example was filtered earlier to have min depth of 8 and now has a max depth of 60)
bcftools view -i 'FORMAT/DP>=8 & FORMAT/DP <60' -O v -o ecology_dp8_sspc_paralog_x0.75_nobadSites_maxdp60.vcf ecology_dp8_sspc_paralog_x0.75_nobadSites.vcf.gz

# Remove replicate samples if wanted to. or any sample you dont want included. -s if you want to list samples. -S if you want to show it a txt file conatining all samples you want exclude. ^ to exclude from new file
bcftools view -s ^B712a,B712b full_snps_dp8_sspc_paralog_x0.75_nobadsites_maxdp60.vcf > full_snps_dp8_sspc_paralog_x0.75_nobadsites_maxdp60_TR.vcf

conda deactivate
logout 



