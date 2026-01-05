#!/bin/bash

# Job specifications for the HPC cluster
#SBATCH --account=def-bureau
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-01:00
#SBATCH --output=/Slurm_files/FounderRare_%A.out
#SBATCH --mail-user=samir.oubninte.1@ulaval.ca
#SBATCH --mail-type=FAIL


  
  
# Example running founderRare approche :

## Prepare phased chip genotype inputs, keep affected probands (IIDs), infer IBD (Identical By Descent) segments
## on the chosen chromosome, remove segments crossing genome gaps, cluster IBD segments, filter long clusters (>10 Mb),
## and compute the Smsg statistic for that chromosome.

# To run code, example :
  #for chrom in 22; do sbatch --export=ALL,chrom=${chrom} run_FounderRare.sh ; done

# Print the current chromosome 
echo "****************************  chr $chrom**************************** " 

# Load necessary modules
module load StdEnv/2020
module load java
module load vcftools
module load gcc/9.3.0
module load r-bundle-bioconductor/3.16
module load r/4.2.1



# Define input paths and create necessary directories
path_work="$(pwd)/Smsg_RealData"
mkdir -p $path_work
path_inputs="$path_work/inputs"
chip_phased="genotypage_GSA_Omni/Union_chr${chrom}_affected.vcf.gz"
affected_list="SCZ.Realaff.HorsFNB.proband.txt"



# Keep only the affected individuals' data from the phased genotype chip data
out_vcf="${path_inputs}/"
mkdir -p $out_vcf
vcftools  --keep $affected_list --gzvcf $chip_phased --recode --out "${out_vcf}/genotype_chip_chr${chrom}.phased.affected"



# Infer identity by descent (IBD) using the refined IBD tool
path_RefinedIBDS="$path_work/RefinedIBD_PIBDS"
mkdir -p "$path_RefinedIBDS"
ibd_file="$path_RefinedIBDS/PSIBD_chr${chrom}.affected"
java -jar /refined-ibd.17Jan20.102.jar gt="${out_vcf}/genotype_chip_chr${chrom}.phased.affected.recode.vcf" length=1.5 out="$ibd_file"



# Delete gaps in the IBD data using a custom R script
gaps_file="gaps_hg19.txt"
ibd_file_in="${ibd_file}.ibd.gz"
ibd_file_out="${ibd_file}.sansGaps.ibd"
Rscript /Psibd_delete_Gaps.r $gaps_file $ibd_file_in $ibd_file_out



# Infer clusters of IBD using Dash-adv, include density and minimum haplotype length
Dash_clusters="$path_work/Dash_clusters/"
mkdir -p "$Dash_clusters"
famfile="genotype_chip_chr${chrom}.fam" 
ibd_file=$ibd_file_out
clusterIBD="$Dash_clusters/IBDCls_chr${chrom}"
Density=0.6
minHaplo=2
cat ${ibd_file} | awk '$9>1.5'| cut -f1,2,3,4,6,7 | sed s/_/' '/g | awk -F " " '{print $1,$2,$3,$4,$5,$6,$7,$8}' |sed s/' 0 '/'.0 '/g | sed s/' 1 '/'.1 '/g | sed s/' 2 '/'.2 '/g | /dash-1-1/bin/dash_adv -fam ${famfile}.fam -density $Density -min $minHaplo "${clusterIBD}_D${Density}_Min${minHaplo}hap"




# Delete clusters longer than 10 Mbp
clusterIBD="${clusterIBD}_D${Density}_Min${minHaplo}hap"
cp ${clusterIBD}.hcl ${clusterIBD}.avecLC.hcl
awk '($3-$2)<10000000' ${clusterIBD}.avecLC.hcl > ${clusterIBD}.hcl
rm ${clusterIBD}.avecLC.hcl




# Calculate Smsg
Smsg="$path_work/Smsg"
mkdir -p "$Smsg"
Rscript Smsg.R clsFileNameWithPath=${clusterIBD}.hcl chr=$chrom  Naff=$(wc -l <${affected_list}) outFileNameWithPath="${Smsg}/Smsg_chr${chrom}"

echo "**************************** ... End **************************** "
