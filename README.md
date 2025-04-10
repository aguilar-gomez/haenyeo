# haenyeo
## Publication associated
Aguilar-Gómez, D., Bejder, J., Andreasen, J., Tristani-Firouzi, M., Ko, Y., Nordsborg, N., Lee, J. Y., Nielsen, R., & Ilardo, M. (In press, Cell Reports.). Adaptation to diving in the Haenyeo breath-hold divers of Korea
https://doi.org/10.1016/j.celrep.2025.115577

## Analysis
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
  
 **5. Structure and Selection**
  - Convert from plink to OHANA format
  - Structure analysis using *qpas*
  - Selection scan using *selscan*
  - Plotting of scan and selection of SNPs
 
**6. Genetic association**
  - Test particular SNPs using *PLINK*

**7. Divergence time and Population Size**
  - Use *RELATE* that uses Ancestral Recombination Graphs
  - Use  *SMC++* that combines the computational efficiency of the SFS and the advantage of utilizing LD information in coalescent HMMs

**8. Testing selection over time**
  - Use *CLUES2* to test selection over time in our candidate SNP

**9. Demographic inference**
  - Use *Treemix* to infer population split and admixture events
  - Use *Admixture Bayes*
