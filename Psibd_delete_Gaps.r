# Load necessary library
library(GenomicRanges)

# Retrieve the command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Get the file paths for the gaps file, input IBD file, and output IBD file
gaps_file=args[1]
ibd_file_in=args[2]
ibd_file_out=args[3]

# Read the gaps file
gaps=read.table(gaps_file, header=T)

# Read the input IBD file
PSIBD.chrom=read.table(ibd_file_in, header=F)
# Get the unique chromosome number from the IBD file
chr=unique(PSIBD.chrom$V5)

# Create a GRanges object for the gaps
gr.gaps=GRanges(
  seqnames = chr ,
  ranges = IRanges(start=gaps$chromStart[gaps$chrom==paste0("chr",chr)],
                   end = gaps$chromEnd[gaps$chrom==paste0("chr",chr)]),
  score = "*")

# Create a GRanges object for the IBD segments
gr.PSIBD.chrom=GRanges(
  seqnames = chr ,
  ranges = IRanges(start = PSIBD.chrom$V6,
                   end = PSIBD.chrom$V7),
  score = "*")

# Find overlaps between the gaps and the IBD segments
ov= findOverlaps(gr.gaps , gr.PSIBD.chrom )

# If there are overlaps, remove the overlapping IBD segments
if (length(ov)!=0) {
  PSIBD.chrom=PSIBD.chrom[-c(unique(subjectHits(ov))),]
} else  {cat("chr :", chr, "length ov :", length(ov))}

# Write the remaining IBD segments to the output file
write.table(PSIBD.chrom, file = ibd_file_out, sep = "\t",
            row.names = FALSE, col.names = F, quote=F)
