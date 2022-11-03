# haenyeo

**1. Preprocessing Reads**
  - Clean reads using *prinseq-lite*
  - Mapping using *bwa-mem*

**2. Imputation**
  - Imputation using 400 individuals from Korea1K Project
  - Using *GLIMPSE*

**3. Merging and filtering**
  - Using *bcftools* and *plink*
  - filter SNPS
  - remove related individuals
  - PCA plot
  
 **4. OHANA**
  - Convert from plink to OHANA format
  - Structure analysis using *qpas*
  - Selection scan using *selscan*
