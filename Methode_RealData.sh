#!/bin/bash

#SBATCH --account=def-girardsi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-00:30

#SBATCH --output=~/projects/rrg-girardsi/genealogy_sims/results/Samir/Real_Vs_Sims/New_RealData2022/Smsg_RealData/Slurm_%A.out
#SBATCH --mail-user=samir.oubninte.1@ulaval.ca
#SBATCH --mail-type=FAIL

# to run cmd example :
  #for chrom in 22; do sbatch --export=ALL,chrom=${chrom} Methode_RealData.sh ; done
  
  
# Example running founderRare approche :

# Prepare phased chip genotype inputs, keep affected probands (IIDs), infer IBD (Identical By Descent) segments
# on the chosen chromosome, remove segments crossing genome gaps, cluster IBD segments, filter long clusters (>10 Mb),
# and compute the Smsg statistic for that chromosome.


echo "****************************  chr $chrom **************************** " 

echo "**************************** Import module ****************************"
module load StdEnv/2020
module load java
module load vcftools
module load gcc/9.3.0
module load r-bundle-bioconductor/3.16
module load r/4.2.1
module load bcftools
module load plink/1.9b_6.21-x86_64



echo "**************************** Inputs and filters ****************************"
path_work="$(pwd)/Smsg_RealData"
mkdir -p $path_work

path_inputs="$path_work/inputs"
mkdir -p $path_inputs

#Files are already phased, filtered by freq (MAF)
chip_phased="${path_inputs}/puces_chr${chrom}_filtred_probands.recode.vcf"

#Check if NB family is exluded (Not yet)
#Give up using fam file 


echo "****************************   Keep affected IIDs **************************** " 
affected_list="$(pwd)/probands_scz_affected2024.txt"
#vcftools  --keep $affected_list --gzvcf $chip_phased --recode --out "${path_inputs}/genotype_chip_chr${chrom}.phased.affected"
vcftools --vcf  $chip_phased --keep $affected_list --recode --out ${path_inputs}/puces_chr${chrom}_filtred_affectedPro



echo "**************************** Infer ibd by refinedIBD  **************************** "
path_RefinedIBDS="$path_work/RefinedIBD_PIBDS"
mkdir -p "$path_RefinedIBDS"
    
ibd_file="$path_RefinedIBDS/PSIBD_chr${chrom}.affectedPro"
java -jar /lustre03/project/6033529/genealogy_sims/results/Samir/software/refined-ibd.17Jan20.102.jar gt="${path_inputs}/puces_chr${chrom}_filtred_affectedPro.recode.vcf" length=1.5 out="$ibd_file"

echo "**************************** Deleting gaps **************************** "
gaps_file="/lustre03/project/6033529/genealogy_sims/results/Samir/software/gaps_hg38.txt"
ibd_file_in="${ibd_file}.ibd.gz"
ibd_file_out="${ibd_file}.sansGaps.ibd"
Rscript /lustre03/project/6033529/genealogy_sims/results/Samir/software/Psibd_delete_Gaps.r $gaps_file $ibd_file_in $ibd_file_out



echo "**************************** Infer Cluster of IBD by Dash-adv **************************** "
Dash_clusters="$path_work/Dash_clusters/"
mkdir -p "$Dash_clusters"
# famfile="/lustre03/project/6033529/genealogy_sims/results/Samir/Sims_power/WGS_to_Chip/Phased_Chip/genotype_chip_chr${chrom}rep1.Good.fam" # *chr${chrom}rep${rep}.fam all this files are same
ibd_file=$ibd_file_out
clusterIBD="$Dash_clusters/IBDCls_chr${chrom}"

# Parameters : 
Density=0.6
minHaplo=2
#PSIBD > 1.5 and input format for Dash
#cat ${ibd_file} | awk '$9>1.5'| cut -f1,2,3,4,6,7 | sed s/_/' '/g | awk -F " " '{print $1,$2,$3,$4,$5,$6,$7,$8}' |sed s/' 0 '/'.0 '/g | sed s/' 1 '/'.1 '/g | sed s/' 2 '/'.2 '/g | /lustre03/project/6033529/SOFT/dash-1-1/bin/dash_adv -fam ${famfile}.fam -density $Density -min $minHaplo "${clusterIBD}_D${Density}_Min${minHaplo}hap"
cat ${ibd_file} | awk '$9>1.5'| cut -f1,2,3,4,6,7 | sed s/_/' '/g | awk -F " " '{print $1,$2,$3,$4,$5,$6,$7,$8}' |sed s/' 0 '/'.0 '/g | sed s/' 1 '/'.1 '/g | sed s/' 2 '/'.2 '/g | /lustre03/project/6033529/SOFT/dash-1-1/bin/dash_adv -density $Density -min $minHaplo "${clusterIBD}_D${Density}_Min${minHaplo}hap"



echo "**************************** Deleting clusters >10Mbp **************************** "
clusterIBD="${clusterIBD}_D${Density}_Min${minHaplo}hap"
cp ${clusterIBD}.hcl ${clusterIBD}.avecLC.hcl
awk '($3-$2)<10000000' ${clusterIBD}.avecLC.hcl > ${clusterIBD}.hcl
rm ${clusterIBD}.avecLC.hcl


echo "**************************** Calculating Smsg **************************** "
#Creat repertory for S_MSG
Smsg="$path_work/Smsg"
mkdir -p "$Smsg"

#Calculating S_mSG
 Rscript /lustre03/project/6033529/genealogy_sims/results/Samir/software/Smsg.R clsFileNameWithPath=${clusterIBD}.hcl chr=$chrom  Naff=$(wc -l <${affected_list}) outFileNameWithPath="${Smsg}/Smsg_chr${chrom}"

# echo "**************************** ... End **************************** "


