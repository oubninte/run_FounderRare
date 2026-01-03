
# This script processes simulated genomic data to identify rare variants in CRG

import numpy as np
import tskit as ts
import pandas as pd
import sys
import os

# ==========================================
# 1. Input Arguments and File Loading
# ==========================================

# Input arguments provided via the command line:
# 1. path_wgs: Path to the directory containing simulation files
# 2. chrom: Chromosome number (integer)
# 3. rep: Simulation replicate number (integer)
path_wgs = sys.argv[1]
chrom = int(sys.argv[2])
rep = int(sys.argv[3])

# Construct the file path for the tree sequence (ts) file
ts_file = (
    f"{path_wgs}genCRG_AndNosF_ForSims.txt.generations.4sims_"
    f"{chrom}_1.66e-8_OutOfAfrica_2T12_rescale_mu1.66_grch38_0.0_sim{rep}.ts"
)

# Load the tree sequence file
try:
    ts = ts.load(ts_file)
except FileNotFoundError:
    print(f"[ERROR] Tree sequence file not found: {ts_file}")
    sys.exit(1)

# ==========================================
# 2. Load CRG and SCZ Individual Data
# ==========================================

# Load the CRG individuals file
gen = pd.read_csv("/lustre03/project/6033529/genealogy_sims/results/Samir/Genealogy/gen.CRG.aveCodesGenetiques_v3_avecGen2.csv")
CRG_ind = gen[pd.notnull(gen["CodeCartagene"])]["ind"]  # Filter CRG individuals by non-null 'CodeCartagene'

# Load the probands file
probands = pd.read_csv("/lustre03/project/6033529/genealogy_sims/results/Samir/Genealogy/Probands_CRGAndNosF_ForSims.csv", header=None)

# Identify CRG individuals present in probands
CRG_ind = np.intersect1d(CRG_ind, probands)
CRG_ind = set(CRG_ind)

# Identify SCZ individuals by subtracting CRG from probands
probands = set(probands[0])  # Convert probands DataFrame column to set
SCZ_ind = probands - CRG_ind

# Convert identifiers to strings for compatibility
CRG_ind = [str(x) for x in CRG_ind]
SCZ_ind = [str(x) for x in SCZ_ind]

# ==========================================
# 3. Simplify Tree Sequence for CRG and SCZ
# ==========================================

# Extract nodes belonging to CRG individuals
Nodes_CRG = []
for ind in ts.individuals():
    if ind.metadata["individual_name"] in CRG_ind:
        Nodes_CRG.append(ind.nodes[0])
        Nodes_CRG.append(ind.nodes[1])

# Simplify the tree sequence for CRG individuals
CRG_ts = ts.simplify(Nodes_CRG)

# Extract nodes belonging to SCZ individuals
Nodes_SCZ = []
for ind in ts.individuals():
    if ind.metadata["individual_name"] in SCZ_ind:
        Nodes_SCZ.append(ind.nodes[0])
        Nodes_SCZ.append(ind.nodes[1])

# Simplify the tree sequence for SCZ individuals
SCZ_ts = ts.simplify(Nodes_SCZ)

# ==========================================
# 4. Split Genome into Bins
# ==========================================

# Define the length of each bin (400,000 base pairs)
length_region = 400_000

# Calculate the number of bins (M)
M = int(SCZ_ts.sequence_length // length_region)
print(f"[INFO] Split genome into {M} bins of length 400K")

# Generate the start positions for each bin
# (bins are indexed from 0 to M-1 because of Python indexing)
l_pot_causal_start = range(0, (M - 1) * length_region, length_region)

# ==========================================
# 5. Identify Rare Variants in Regions
# ==========================================

# Lists to store results
N_VR_CRG = []  # Number of rare variants in CRG
N_VnonR_SCZ = []  # Number of non-rare variants in SCZ
N_V_SCZ = []  # Number of variants in SCZ
start_pos_l = []  # Start positions of bins

for start_pos in l_pot_causal_start:
    end_pos = start_pos + length_region

    # Restrict the CRG tree sequence to the current region
    restricted_CRG_ts = CRG_ts.keep_intervals([[start_pos, end_pos]])
    MAF = []  # Minor allele frequencies
    Position = []  # Positions of variants
    RA = []  # Rare alleles

    # Process variants in the restricted CRG tree sequence
    for v in restricted_CRG_ts.variants():
        site = v.site
        alleles = np.array(v.alleles)
        values, counts = np.unique(alleles[v.genotypes], return_counts=True)

        # Skip tri-allelic sites
        if len(values) > 2:
            continue

        # Identify rare variants (MAF < 0.01)
        maf = round(min(counts) / sum(counts), 5)
        if maf < 0.01:
            MAF.append(maf)
            Position.append(site.position)
            RA.append(values[np.where(counts == min(counts))][0])

    # Store rare variants for CRG
    RARE_VARIANT = {"Position": Position, "MAF": MAF, "RA": RA}

    # Restrict the SCZ tree sequence to the current region
    SCZ_ts_restricte = SCZ_ts.keep_intervals([[start_pos, end_pos]])

    # Process variants in the restricted SCZ tree sequence
    MAF_SCZ = []
    Position_SCZ = []
    RA_SCZ = []
    for v in SCZ_ts_restricte.variants():
        site = v.site

        # Only consider positions present in CRG rare variants
        if site.position not in RARE_VARIANT["Position"]:
            continue

        alleles = np.array(v.alleles)
        values, counts = np.unique(alleles[v.genotypes], return_counts=True)
        f = round(min(counts) / sum(counts), 5)

        # Skip variants with frequency > 0.5
        if f > 0.5:
            continue

        MAF_SCZ.append(f)
        Position_SCZ.append(site.position)
        RA_SCZ.append(values[np.where(counts == min(counts))][0])

    # Store rare variants for SCZ
    VR_MAF_inSCZ = {"Position": Position_SCZ, "MAF": MAF_SCZ, "RA": RA_SCZ}

    # Count variants based on MAF thresholds
    N_sup = sum(1 for maf in VR_MAF_inSCZ["MAF"] if maf >= 0.01)
    N_sup2 = sum(1 for maf in VR_MAF_inSCZ["MAF"] if maf >= 0.03)

    # Append results for the current bin
    N_VR_CRG.append(len(RARE_VARIANT["Position"]))
    N_V_SCZ.append(len(VR_MAF_inSCZ["Position"]))
    N_VnonR_SCZ.append(N_sup2)
    start_pos_l.append(start_pos)

# ==========================================
# 6. Save Results to CSV
# ==========================================

# Create a DataFrame to store the results
df0 = pd.DataFrame({
    "Start": start_pos_l,
    "N_VR_CRG": N_VR_CRG,
    "N_V_SCZ": N_V_SCZ,
    "N_VnonR_SCZ": N_VnonR_SCZ
})

# Define the output path
output_dir = "/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/Sims_power/Phenotype_sims_100reps/ForAlexandre/"
output_file = f"{output_dir}/selection_region_chr{chrom}rep{rep}.csv"

# Save the DataFrame to a CSV file
df0.to_csv(output_file, index=False)
#print(f"[INFO] Results saved to {output_file}")
