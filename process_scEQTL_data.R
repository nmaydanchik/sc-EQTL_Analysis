# returns a data frame of cell-type eQTL and Alzheimer's Disease GWAS, already LD Clumped
# Needs Y (Alzheimer's Disease GWAS) and SNP_Alleles to be loaded.

process_sc_eQTL_Data<- function(filenum, study) {
  if (study == "Bryois") filename = paste0("C:/Users/Nate/Documents/RLibrary/Bryois/Bryois2022NN_",filenum,".sig_qtl.tsv")
  if (study == "Jerber") filename = paste0("C:/Users/Nate/Documents/RLibrary/Jerber/Jerber2021NG_",filenum,".sig_qtl.tsv")
  
  X <- fread(filename, select = c("geneName", "variantId", "pValue", "beta", "se", "cellTypeName")) # Load eQTL data
  
  # Filter for the CAV1 gene
  X <- X %>% 
    filter(geneName == "CAV1")
  colnames(X) = c("geneName", "rsid", "pval", "beta", "se", "cellClusterName")
  
  if (nrow(X) > 0) {
    
    # Attach alleles to eQTL data
    X <- inner_join(X, SNP_Alleles, by = join_by(rsid == SNP)) %>% 
      rename("non_effect_allele" = "other_allele")
    
    # Merge with AD data
    df <- inner_join(X, Y, by = join_by(rsid))
    
    # Match alleles between the eQTL and GWAS data
    for (i in 1:nrow(df)) {
      if (df[i,]$effect_allele == df[i,]$non_effect_allele.ad && df[i,]$non_effect_allele == df[i,]$effect_allele.ad) {
        df[i,]$beta = -1 * df[i,]$beta # If they are switched, invert beta of eQTL
      }
      else if (df[i,]$effect_allele != df[i,]$effect_allele.ad || df[i,]$non_effect_allele != df[i,]$non_effect_allele.ad)  {
        print(paste0("Allele mismatch: file ", filenum, ", row ", i)) # If they are not switched and not the same, output a warning
      }
    }
    
    # filter by p-value and LD clump
    df <- df %>% select(!(c("effect_allele.ad", "non_effect_allele.ad", "effect_allele", "non_effect_allele"))) %>% 
      filter(pval < 0.005) %>% 
      ld_clump(clump_kb = 250, clump_r2 = 0.5) %>% 
      select(!"id")
  }
} 
