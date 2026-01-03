

#  Objective:
## Simulate case/control phenotypes per chromosome and replicate for each rare variant in specified regions.
## Use a logistic model with baseline prevalence 1% and apply an odds ratio (OR = `odd`) to one causal variant at a time.
## Count affected individuals per family (from `pedsimple2`) and select variants that meet recruitment rule:
## at least 2 affected per family and at least 16 affected overall; save the affected IDs for the best variant per region.
##


rm(list=ls())
library(stringr)

### Import data#######
# chr=21
# rep=12
# odd=30

args <- commandArgs(trailingOnly = TRUE)
chr = as.numeric(args[1])
rep = as.numeric(args[2])
odd = as.numeric(args[3])
chemin_save = args[4]


chemin.VR="/lustre03/project/6033529/genealogy_sims/results/Samir/Sims_power/simulate_pheno/Rare_Variant_sims/"
regions=read.table(paste0("/lustre03/project/6033529/genealogy_sims/results/Samir/Sims_power/simulate_pheno/", "ID_regions.txt"), header = T, sep="\t")
regions=regions[regions$chr==chr &regions$rep==rep,]
regions=regions[!duplicated(regions),]

for (i in regions$Index) {
  
  RareVariant= read.table(paste0(chemin.VR, "Rare_Variant_chr",chr,"rep",rep,"_region",i,".csv") , header = T, sep=",")
  #set.seed(125)
  Matrice_Genotype=RareVariant[,-c(1)] 
  N_vr=ncol(Matrice_Genotype)
  
  #### simulation de phenotype : modele de base #########
  
  #-------------------------------------------------------------
  #### une seule variante causal #########
  OR =c(1,odd)
  beta0=log(0.01/(1-0.01))
  
  set.seed(123450)
  #Selection de VR
  xx=apply(Matrice_Genotype, 2, function(l){
    sum(l>0)
  }  )
  cat("\n la variant qui la plus portée", names(xx[ which.max(xx)]) ,"\n nombre de porteur ", xx[ which.max(xx)], "\n chr",chr,"rep",rep,"region",i,"odd",odd,"\n")
  #causalVariantIndex= which.max(xx)# sample(1:ncol(Matrice_Genotype),1) #10 #which.max(xx)
  
  Recap_allVR=data.frame()
  pheno_allVR=data.frame()
  load("/lustre03/project/6033529/genealogy_sims/results/Samir/genotypage_GSA_Omni/RVS/pedsimple2.RData")
  
  #do for every rare variant
  for (j in 1:ncol(Matrice_Genotype)) {
    causalVariantIndex=j
    pheno=data.frame(RareVariant[,1])
  
    #logistic model
    beta=rep(log(OR[1]), N_vr)
    beta[causalVariantIndex]=rep(log(OR[2]), length(causalVariantIndex))
    pY=1/(1+exp(-beta0- as.matrix(Matrice_Genotype) %*% beta))
    Y=rbinom(length(pY) , 1, pY) 
    pheno=cbind(pheno, round(Y))
    
    colnames(pheno)=c("ID", "Y")
    Recap=data.frame(OR=odd)
    
    # import families ..
    famid=unlist(sapply(pheno[,1], function(ind)    pedsimple2$famid[pedsimple2$id %in% ind]))
    pheno=pheno[pheno$ID %in% pedsimple2$id, ]
    pheno=cbind(pheno, data.frame(famid))
    
    # nombre de sujets affectés
    Recap$N_Aff=sum(pheno$Y)
    
    # nombre familles contenant au moins 1 affecté
    Recap$"NFam_1+"=length(unique(pheno$famid[pheno$Y==1]))
    
    # nombre familles contenant au moins 2 affec
    Recap$"NFam_2+"=sum(table(pheno$famid[pheno$Y==1])>=2)
    
    # Save infos about VR
    Recap$VR=(colnames(Matrice_Genotype)[j])
    Recap$IndexVR=j
    Recap_allVR=rbind(Recap_allVR,Recap)
    
    pheno$VR=(colnames(Matrice_Genotype)[j])
    pheno$IndexVR=j
    pheno$OR=OR[2]
    pheno_allVR=rbind(pheno_allVR,pheno)
  }
  
  #save affected list that respect condition min 16 affected and min 2 affected per family
  Just_affected=pheno_allVR[pheno_allVR$Y==1 & pheno_allVR$OR==odd,]
  Just_affected2=data.frame()
  
  # min 2 affected per family
  for (VR in unique(Just_affected$VR)) {
    keep=table(Just_affected$famid[Just_affected$VR == VR])>=2
    fam2keep=names(table(Just_affected$famid[Just_affected$VR == VR]))[keep]
    Just_affected2=rbind(Just_affected2,Just_affected[Just_affected$famid %in% fam2keep & Just_affected$VR == VR,])
  }

    # min 16 affected
  VR2keep=names(table(Just_affected2$VR))[table(Just_affected2$VR)>=16]
  # save VR that give max of affected 
  if (length(VR2keep)!=0) {
    Vrt= names(which(table(Just_affected2$VR)==max(table(Just_affected2$VR))))[1] # variant qui donne le max d'affectés
    affecteds=paste0(Just_affected2$famid[Just_affected2$VR==Vrt],"_",Just_affected2$ID[Just_affected2$VR==Vrt])
  }

write.table( affecteds ,paste0(chemin_save,"affectedList_chr",chr,"rep",rep,"_region",i,"OR",odd,".txt"),sep="\t", quote=F, row.names = F, col.names = F)
}
