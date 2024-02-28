# Clear all variables from the workspace
rm(list = ls())

# Load necessary libraries
library(FounderRare)
library(Rcpp)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(igraph)
library(stringr)

# Function to process command line arguments
process_args <- function(args_input) {
  # Split the input arguments by space
  args1=strsplit(args_input, " ")
  # Split the arguments by "=" to separate names and values
  args2=str_split_fixed(args1,"=",length(args1))
  # Extract the argument values
  args=args2[,2]
  # Assign names to the argument values
  names(args)=args2[,1]
  return(args)
}

# Retrieve arguments
args_input <- commandArgs(trailingOnly = TRUE)
args=process_args(args_input)
chr  = as.numeric(args["chr"])
clsFileNameWithPath=args["clsFileNameWithPath"]

# Check if Naff is provided in the arguments
if ("Naff" %in% names(args)) {
  # Convert Naff to numeric
  Naff=as.numeric(args["Naff"])
  # Read the cls file with column names based on Naff
  IBDClusters=read.table(clsFileNameWithPath, col.names = paste0("V",seq_len(2*2*Naff+3)), fill=T)
} else{
  # If Naff is not provided, assume a default value of 10000 (the number of people affected must be less than 1000)
  IBDClusters=read.table(clsFileNameWithPath, col.names = paste0("V",seq_len(2*2*10000+3)), fill=T)
}

# Check if output file path is provided in the arguments
if ("outFileNameWithPath"  %in% names(args) ) {
  outFileNameWithPath = args["outFileNameWithPath"]
} else {
  # If not provided, use a default output file name
  outFileNameWithPath="Results"
}

# Start inferring SG
cat("\n Inferring SG ...  \n")
# Remove columns with all NA values
IBDClusters=IBDClusters[, colSums(is.na(IBDClusters)) != nrow(IBDClusters)]
# Call the TBM function to infer SG
out=TBM(IBDClusters, chr=chr)
# Save the output to a RData file
saveRDS(out, file=paste0(outFileNameWithPath, "_SGList.RData"))

# Start inferring SG with disjoint clusters
cat("\n SG with disjoint clusters ...  \n")
# Call the TBMD function to infer SG with disjoint clusters
outWithDisjClusters=TBMD(out, IBDClusters)
# Save the output to a RData file
saveRDS(outWithDisjClusters, file=paste0(outFileNameWithPath, "_SGList_CD.RData"))

# Start computing Smsg
cat("\n Computing Smsg ...  \n")
# Check if Naff is provided in the arguments
if ("Naff" %in% names(args)) {
  # If Naff is provided, call the Smsg function with Naff
  outSmsg=Smsg(outWithDisjClusters, Naff)
} else {
  # If Naff is not provided, call the Smsg function without Naff
  outSmsg=Smsg(outWithDisjClusters)
}

# Write the output to a csv file
write.table( outSmsg , paste0(outFileNameWithPath, ".csv") ,sep=",", quote=F, row.names = F, col.names = T)
