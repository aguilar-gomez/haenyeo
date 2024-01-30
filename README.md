# haenyeo

**1. Preprocessing Reads**
  - Clean reads using *prinseq-lite*
  - Mapping using *bwa-mem*

**2. Imputation and Phasing**
  - Imputation using 400 individuals from Korea1K Project
  - Using *GLIMPSE*

**3. Relatedness**
  - Calculate relatedness
  - Remove related individuals

**4. Merging and filtering**
  - Using *bcftools* and *plink*
  - filter SNPS
  - remove related individuals
  - maximum amount of genotypes missing 
  - remove positions with change in reference allele
  - PCA plot
  
 **5. OHANA**
  - Convert from plink to OHANA format
  - Structure analysis using *qpas*
  - Selection scan using *selscan*

 
**6. Plotting and selecting SNPs for association**
  - Plotting of scan and selection of SNPs
 
**7. Genetic association**
  - Test particular SNPs using *gemma*

**8. Divergence time and Population Size**
  - Use *RELATE* that uses Ancestral Recombination Graphs
  - Use  *SMC++* that combines the computational efficiency of the SFS and the advantage of utilizing LD information in coalescent HMMs
