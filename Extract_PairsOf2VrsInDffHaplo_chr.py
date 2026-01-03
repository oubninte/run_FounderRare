#-------------------------------------------------------------------------------
# Name:        VR in diff haplotype
# Purpose:
#
# Author:      SO
#
# Created:     05-09-2023
# Copyright:   (c) Saoub9 2023
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#In this code, I store the pairwise VR in different haplo,
#for example, in the output file, the first and the second Vr are in different haplo
#The 3e and the 4e are in different haplo, but not sur if the 2e and the 3e are in diff haplo

def main():
    pass

if __name__ == '__main__':
    main()


import numpy as np
import tskit as ts
import pandas as pd
import sys
import os


#***********************************************************************************************************************************
# 1. input argument for this script
#***********************************************************************************************************************************

chrom= int(sys.argv[1])
rep= int(sys.argv[2])
N_region= int(sys.argv[3])
chemin_RV=sys.argv[4]
chemin_tsFile=sys.argv[5]
chemin4save=sys.argv[6]

#For test
#chrom=21
#rep= 100
#N_region=1
#chemin_RV="/lustre03/project/6033529/genealogy_sims/results/Samir/Sims_power/Phenotype_sims_100reps/Rare_variants/"
#chemin_tsFile="/lustre03/project/6033529/genealogy_sims/results/Samir/Sims_data/"
#chemin4save=""

regions=pd.read_csv("/lustre03/project/6033529/genealogy_sims/results/Samir/Sims_power/Phenotype_sims_100reps/regions_casaux_GOOD.txt", delimiter=" ")
regions=regions.drop_duplicates()

#***********************************************************************************************************************************
# 2. Separate CaG and SCZ and get 
#***********************************************************************************************************************************

#Vrs file and there genotype
Rare_Variant=pd.read_csv(chemin_RV+"Rare_Variant_chr"+str(chrom)+"rep"+str(rep)+"_region"+str(N_region)+".csv", delimiter=",")
#np.sum(Rare_Variant["VR322"])

#Extract variant by there site.id in column names of Rare_Variant df (take num-1)
ts_file=chemin_tsFile+"genCRG_AndNosF_ForSims.txt.generations.4sims_"+str(chrom)+"_1.66e-8_OutOfAfrica_2T12_rescale_mu1.66_grch38_0.0_sim"+str(rep)+".ts"
ts=ts.load(ts_file)

#TS of CacG and TS of SCZ
    # firting SCZ and CRG individuals
chemin_gen="/lustre03/project/6033529/genealogy_sims/results/Samir/Genealogy/"
gen=pd.read_csv(chemin_gen+"gen.CRG.aveCodesGenetiques_v3_avecGen2.csv")
CRG_ind=gen[ pd.notnull(gen["CodeCartagene"])]["ind"]
probands=pd.read_csv(chemin_gen+"/Probands_CRGAndNosF_ForSims.csv", header=None)

CRG_ind=np.intersect1d(CRG_ind,probands)
CRG_ind=set(CRG_ind)

#probands
probands= set(probands[0])# probands[0].values.tolist()
SCZ_ind=probands-CRG_ind

# Sans la converstion, le ts sera vide
SCZ_ind=[str(x) for x in SCZ_ind]


# Se reistreindre au sujets SCZ pour cacluler MAF
Nodes_SCZ=[]
for ind in ts.individuals():
    if ind.metadata["individual_name"] in SCZ_ind :
        Nodes_SCZ.append(ind.nodes[0])
        Nodes_SCZ.append(ind.nodes[1])

SCZ_ts = ts.simplify(Nodes_SCZ)

start_pos=int(regions[(regions["chr"]==chrom) & (regions["rep"]==rep) & (regions["Index"]==N_region)]["Start"])
end_pos=int(regions[(regions["chr"]==chrom) & (regions["rep"]==rep) & (regions["Index"]==N_region)]["End"])
#ID of individuals (to keep the same order in ts other way i have already SCZ_ind, the order is important)
SCZ_ts_restricte = SCZ_ts.keep_intervals([[start_pos,end_pos]])
sample_nodes = SCZ_ts_restricte.samples()
sample_individuals = []
for ind in SCZ_ts_restricte.individuals():
    if len(ind.nodes) == 0:
        continue
    # diploid - two nodes, both the same individual
    if ind.nodes[0] in sample_nodes:
        assert len(ind.nodes) == 2
        assert ind.nodes[1] in sample_nodes
        sample_individuals.append(ind)

ID = [  str(ind.metadata["individual_name"]) for ind in sample_individuals ]

column_names = Rare_Variant.columns
VR_siteID=[(int(name[2:])-1) for name in column_names[1:]]  #first col is not a VR


#Get the site.id of VR from the file site.id != site.position
column_names = Rare_Variant.columns
VR_siteID=[(int(name[2:])-1) for name in column_names[1:]]  #first col is not a VR

#Get genotype of these Variant
df=pd.DataFrame({'ID':ID})
for v in SCZ_ts_restricte.variants():
    site = v.site

    if site.id not in VR_siteID:
        continue

    alleles = np.array(v.alleles)
    values, counts = np.unique(alleles[v.genotypes],return_counts=True)

    G=[]
    for r in range(0, len(v.genotypes),2):
        G.append(v.genotypes[r]+v.genotypes[r+1])
    df = pd.concat([df, pd.DataFrame({'VR'+str(site.id+1):G})], axis=1)
    #df.to_csv(chemin_data+"Rare_Variant_chr"+str(chrom)+"rep"+str(rep)+"_region"+str(N_region)+"latest.csv", index=False)

#***********************************************************************************************************************************
# 3. Keep just the VR in different haplotype if its exist
#***********************************************************************************************************************************

i=1
VR_siteID_inDiffHaplo=[]
for v in SCZ_ts_restricte.variants():
    site = v.site
    if site.id not in VR_siteID[:-1]:
        continue

    for v2 in SCZ_ts_restricte.variants():
        site2 = v2.site
        if site2.id not in VR_siteID[i:]:
            continue

        deleteFirstOne=[]
        for r in range(0, len(v.genotypes),2):
            #True if in same haplotype for this individual, otherwise False
            deleteFirstOne.append((v.genotypes[r] !=v2.genotypes[r]) and (v.genotypes[r]!=v.genotypes[r+1]) and (v.genotypes[r] ==v2.genotypes[r+1]) )

        if sum(deleteFirstOne)!=0 :
            VR_siteID_inDiffHaplo.extend([site.id, site2.id])

    i=i+1


# 4. saving
if len(VR_siteID_inDiffHaplo)!=0:
  df=pd.DataFrame({'ID':ID})
  for VrSite in VR_siteID_inDiffHaplo:
      for v in SCZ_ts_restricte.variants():
          site=v.site
          if site.id !=VrSite :
              continue
  
          alleles = np.array(v.alleles)
          values, counts = np.unique(alleles[v.genotypes],return_counts=True)
  
          if min(counts) != sum(v.genotypes):
              print("delete ", site.position ,values,counts, min(counts), sum(v.genotypes))
              continue
  
          G=[]
          for r in range(0, len(v.genotypes),2):
              G.append(v.genotypes[r]+v.genotypes[r+1])
          df = pd.concat([df, pd.DataFrame({'VR'+str(site.id+1):G})], axis=1)
      
  df.to_csv(chemin4save+"PairsOf2VrsInDffHaplo_chr"+str(chrom)+"rep"+str(rep)+"_region"+str(N_region)+".csv", index=False)
  
else:
    print("No couple of variants in different haplotypes.")
  
  
