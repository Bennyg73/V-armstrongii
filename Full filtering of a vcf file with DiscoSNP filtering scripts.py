#full filtering of a clustered vcf file. 
# first start with depth, missing data and Minor allele frequency 
 python3 filter_vcf_by_indiv_cov_max_missing_and_maf.py -i discoRad_k_31_c_3_D_0_P_5_m_5_clustered.vcf -o ecology_c_8_f_0.05_m_0.vcf -c 8 -f 0.05 -m 0 -s 1
#     OUTPUT
# filter_vcf_by_indiv_cov_max_missing_and_maf.py
# input_file : discoRad_k_31_c_3_D_0_P_5_m_5_clustered.vcf
# output_file : ecology_c_8_f_0.05_m_0.vcf
# filter parameters : indiv_DP>=8, missing<=0 (prop*n_geno = 0.0 * 83), snp-only=1
# 2545929 seen variants, 24785 variants after filtering
# initial missing count (on snps only if snp-only=1) = 148782220, genotypes changed to NA = 25607726, final missing count in new vcf = 0
# initial missing percent (on snps only if snp-only=1) = 70.41 %, final missing percent in new vcf = 0.00 %

#then single snppercluster script
python3 SSPC_script.py -i ecology_c_8_f_0.05_m_0.vcf -o ecology_dp8_sspc.vcf

# then paralog filter 
python3 filter_paralogs.py.py -i ecology_dp8_sspc.vcf -o ecology_dp8_sspc_paralog_x0.75.vcf -x 0.75 -y 1
#    OUTPUT
11964 on 13011 clusters had less than 100.0% of SNP with less than 75.0% heterygous genotypes

#to be able to check bad sites file needs to be bgzip and then tabix to extract these bad sites
bgzip ecology_dp8_sspc_paralog_x0.75.vcf

tabix -p vcf ecology_dp8_sspc_paralog_x0.75.vcf.gz

#identify bad SNP sites that may be present.
bcftools view -S sample_file.txt ecology_dp8_sspc_paralog_x0.75.vcf.gz |  bcftools query -i 'N_PASS(GT="AA")<3 & N_PASS(GT="RR")<3 & N_PASS(GT="het")<3' -f '%CHROM\t%POS\n' > bad_sites2.tsv
#double check if the script worked plus to see how many bad sites there are in the file.
wc -l bad_sites2.tsv
2506 bad_sites2.tsv

#Excluding bad sites from our now filtered VCF file.
bcftools view -T ^bad_sites2.tsv -O z -o ecology_dp8_sspc_paralog_x0.75_nobadSites.vcf.gz ecology_dp8_sspc_paralog_x0.75.vcf.gz

# Formating filtered vcf file to have a lower and upper limit for depth (this example was filtered earlier to have min of 8 and now has up to 100 change value to desired max depth), | WC -l is used to tell you how many SNPS are left, instead use > output.vcf to get a file 
# To view amount of SNPS in this format
bcftools query -i 'FORMAT/DP>=8 & FORMAT/DP<100' -f '%CHROM\t%POS\t[%SAMPLE:%DP\t]\n' ecology_dp8_sspc_paralog_x0.75_nobadSites.vcf.gz | wc -l
9311
# To print this format to a new VCF file. 
bcftools view -i 'FORMAT/DP>=8 & FORMAT/DP <60' -O v -o ecology_dp8_sspc_paralog_x0.75_nobadSites_maxdp60.vcf ecology_dp8_sspc_paralog_x0.75_nobadSites.vcf.gz

# Remove replicate samples if wanted to. or any sample you dont want included. -s if you want to list samples. -S if you want to show it a txt file conatining all samples you want exclude. ^ to exclude from new file
bcftools view -s ^B712a,B712b full_snps_dp8_sspc_paralog_x0.75_nobadsites_maxdp60.vcf > full_snps_dp8_sspc_paralog_x0.75_nobadsites_maxdp60_TR.vcf