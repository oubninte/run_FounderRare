## run_FounderRare

This repository hosts the **analysis scripts, pipelines, and reproducibility assets** used in the study
 
>**FounderRare** is a statistical framework and R package to detect genomic regions likely to harbor **rare variants** underlying complex diseases by leveraging **population genealogies** with **founder effects** and **identity‑by‑descent (IBD)** segments inferred from **SNP array** data. The approach partitions the genome into **Synthetic Genomic regions (SGs)**, builds **IBD clusters** of affected individuals, and evaluates **IBD‑sharing enrichment** using (**S_msg**, **S_all**) statistics.

> **Main preprint:**  
> *Statistical Approach Leveraging Genealogies of Populations with a Founder Effect and Identical by Descent Segments to Identify Rare Variants in Complex Diseases.* medRxiv (Sep 18, 2025), DOI: [10.1101/2025.09.16.25335588](https://doi.org/10.1101/2025.09.16.25335588)  
> **Package preprint:**  
> *FounderRare: A Novel Statistical Package to Identify Rare Variants in Complex Diseases.* medRxiv (Sep 11, 2025),  [DOI](https://www.medrxiv.org/content/10.1101/2025.09.11.25335516v1)

---

## Why FounderRare approche?

- **Founder‑population focus.** Incorporates deep **genealogies** (e.g., Eastern Québec families) to exploit relatedness typical of founder effects.
- **Array‑based scalability.** enabling large sample sizes.
- **IBD‑enrichment tests.** Introduces **S_msg** (most‑shared haplotype within an SG) and an adapted **S_all** to quantify whether affected individuals share IBD more than expected under the null.
- **Genealogy‑aware calibration.** Obtains the **genome‑wide maximal statistic’s null** via whole‑genome transmission simulations in the pedigree using **msprime**.

See details in the medRxiv preprints:  
- Main method and results: [10.1101/2025.09.16.25335588](https://doi.org/10.1101/2025.09.16.25335588) (PDF: [full text](https://www.medrxiv.org/content/10.1101/2025.09.16.25335588v1.full.pdf))  
- Package overview: [FounderRare preprint](https://www.medrxiv.org/content/10.1101/2025.09.11.25335516v1)

---

- **FounderRare R package (separate repo)**  
  The production implementation of the **S_msg / S_all** statistics, windowing utilities, cluster building lives at:  
  **GitHub:** https://github.com/oubninte/FounderRare  
  **Package paper:** *FounderRare: A Novel Statistical Package to Identify Rare Variants in Complex Diseases* (medRxiv). [DOI](https://www.medrxiv.org/content/10.1101/2025.09.11.25335516v1)
