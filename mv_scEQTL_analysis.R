library(tidyverse)
library(stringr)
library(data.table)
library(ieugwasr)
library(MendelianRandomization)
setwd("C://Users/Nate/Documents/RLibrary")

# gene_set <- c("NOS3","MMP1", "MMP2", "MMP3", "MMP7", "MMP8", "MMP9", "MMP10", "MMP11", "MMP12", "MMP13", "MMP14", "MMP15", "MMP16", "MMP17", "MMP19", "MMP20", "MMP23A", "MMP23B", "MMP24", "MMP25", "MMP26", "MMP28","NOTCH3","ADAM17","COL4A1","COL4A2","TIMP3","HTRA1","TGFB1","TGFB2","TGFB3","LAMA1","LAMA2","LAMA3","LAMA4","LAMA5","DDR1","DDR2","SOX1", "SOX2", "SOX3", "SOX4", "SOX5", "SOX6", "SOX7", "SOX8", "SOX9", "SOX11", "SOX12", "SOX13", "SOX14", "SOX15", "SOX17", "SOX18", "SOX21", "SOX30", "ITGA1", "ITGA8", "ITGA11", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "NOSTRIN", "MPO", "NOS1", "NOS2", "APP", "PSEN1", "PSEN2", "APOE", "ADAM10", "BCHE", "CLU", "SLC2A1", "ABCB1", "SLC7A5", "CAV1", "IL6", "TNF", "TNFRSF1B", "TNFRSF11B", "PF4", "VWF", "F3", "TFPI", "PLAT", "SERPINE1", "F2", "FGA", "FGB", "FGG", "THBD", "EDN1")
gene_lookup <- read_tsv("gene_info_GTEx_v8.txt.zip") %>% select(c(2, 6)) # Reference table for Gene Symbol and Chromosome
ref <- fread("reference.txt") # Binary table with genes to analyze and AD-related pathways
gene_hits <- lapply(ref, function(X) X==1)

# Load Alzheimer's Disease data
Y <- fread("ieu_AD.txt.gz") %>% 
  select(!"chr")
colnames(Y) = c("non_effect_allele.ad", "effect_allele.ad", "beta.ad", "se.ad", "pval.ad", "rsid")

SNP_Alleles <- fread("snp_pos.txt.gz", select = c("SNP", "effect_allele", "other_allele")) # alleles for sc-eQTL data

for (path in 1:7) {  
  
  gene_output <- gene_lookup[unlist(gene_hits[path+1]),] %>% mutate(numExp = NA, numIV = NA)
  gene_output$Chr = sapply(gene_output$Chr, function(X) substring(X, 4))
  
  for (chrom in 1:22) {
    
    print(paste0("Chromosome ", chrom, ":"))
    
    # Load and transform sc-eQTL Data
    # In my computer, I had a folder "sc_eQTL_AllCells" with 8 folders for each cell type, and then each folder had 22 separate files for each chromosome
    filenames = list.files("sc_eQTL_AllCells", full.names = TRUE, recursive = TRUE, pattern = paste0("\\.",chrom,"\\.gz"))
    celltypes = str_match(filenames, ".*/(.*)/.*")[,2]
    QTL_df <- lapply(filenames, fread, select = c(1,2,4,5)) # QTL_df is a list of data frames, one for each cell type
    QTL_df <- lapply(QTL_df, setNames, c("Gene_id", "rsid", "pval", "beta"))
    names(QTL_df) <- celltypes
    QTL_df <- lapply(QTL_df, mutate, se = abs(beta/qnorm(pval/2))) # Add se to data
    for (i in 1:8) {
      colnames(QTL_df[[i]]) = c("Gene_id", "rsid", paste0("pval.", i), paste0("beta.", i), paste0("se.", i))
    }
    
    for (gene in 1:nrow(gene_output)) {
      if (gene_output[gene, 1] == chrom) {
        # Filter for gene and only keep cell types that have SNPs for that gene
        X <- lapply(QTL_df, function(df) df[str_detect(df$Gene_id, paste0("^", toString(gene_output[gene, 2]), "[\\._]"))]) %>% 
          keep(function(df) nrow(df) > 0)
        
        if(length(X) > 1) {
          # Combine X into one data frame, merge with Alzheimer's Disease data
          df <- X %>% reduce(full_join, by = c("rsid", "Gene_id")) %>% 
            inner_join(SNP_Alleles, by = join_by(rsid == SNP)) %>% 
            rename("non_effect_allele" = "other_allele") %>% 
            inner_join(Y, by = join_by(rsid))
          
          # Match alleles of exposures and outcome. All of the exposures have the same alleles in the data I used.
          for (i in 1:nrow(df)) {
            if (df[i,]$effect_allele == df[i,]$non_effect_allele.ad && df[i,]$non_effect_allele == df[i,]$effect_allele.ad) {
              df[i,]$beta.ad = -1 * df[i,]$beta.ad # If they are switched, invert beta of AD
            }
            else if (df[i,]$effect_allele != df[i,]$effect_allele.ad || df[i,]$non_effect_allele != df[i,]$non_effect_allele.ad)  {
              print(paste0("Allele mismatch: row ", i)) # If they are not switched and not the same, output a warning
              # Add a line here that throws out rows with this problem, if necessary. 
            }
          }
          df <- select(df, !effect_allele:effect_allele.ad)
          
          # Create a temporary data frame that only has p-value columns
          pvals <- select(df, colnames(df)[str_detect((colnames(df)), 'pval')]) %>% 
            select(!pval.ad)
          
          # Filter by p-value. Take any rows that have at least 1 p-value less than threshold
          keptRows = apply(pvals, 1, function(pvals) any(pvals<0.005)) # Some cell types have many, a few hundred, some have few, like 2 or 30
          pvals <- pvals[keptRows]
          
          # Filter df by p-value, add a column that takes the minimum p-value and uses that for LD clumping
          df <- df[keptRows] %>% mutate(pval = apply(pvals, 1, min))
          
          success <- tryCatch(
            {
              df <- df %>% ld_clump(clump_kb = 50, clump_r2 = 0.1) %>% select(!id)
              TRUE
            },
            error = function(e) {
              FALSE
            },
            warning = function(w) {
              FALSE
            }
          )
          if (!success) {
            print(paste0(gene_output[gene, 2], ": LD Failed")) # This would happen because I used the API instead of local.
            next
          }
          
          if (nrow(df) > length(X) && sum(sapply(select(df, colnames(df)[str_detect((colnames(df)), 'pval')]), function(X) sum(X < 0.005))[c(1:length(X))] == 0) == 0) { # Checking that at least as many SNPs as exposures, and that no exposures have 0 significant SNPs.
            # Create matrixes to use mr_mvinput
            df_beta <- select(df, colnames(df)[str_detect((colnames(df)), 'beta')]) %>% 
              select(!beta.ad)
            df_se <- select(df, colnames(df)[str_detect((colnames(df)), 'se')]) %>% 
              select(!se.ad)
            
            # Construct output and run methods
            est = matrix(nrow = length(X), ncol = 8)
            colnames(est) = c("cell_type", "numIV", "Estimate_IVW", "Pvalue_IVW", "Estimate_Egger", "Pvalue_Egger", "Estimate_Median", "Pvalue_Median")
            est[,1] = names(X)
            est[,2] = sapply(select(df, colnames(df)[str_detect((colnames(df)), 'pval')]), function(X) sum(X < 0.005))[c(1:length(X))]
            
            mvmr.obj <- mr_mvinput(bx = as.matrix(df_beta), bxse = as.matrix(df_se), by = df$beta.ad, byse = df$se.ad, exposure = names(X), outcome = "AD")
            res <- mr_mvivw(mvmr.obj)
            est[,3] = res$Estimate
            est[,4] = res$Pvalue
            res <- mr_mvegger(mvmr.obj)
            est[,5] = res$Estimate
            est[,6] = res$Pvalue.Est
            res <- mr_mvmedian(mvmr.obj)
            est[,7] = res$Estimate
            est[,8] = res$Pvalue
            
            write_csv(as.data.frame(est), paste0("scEQTL_output_pathed/", colnames(ref)[path+1], "/", gene_output[gene, 2],".csv"))
            gene_output[gene, 3] = length(X)
            gene_output[gene, 4] = nrow(df)
          }
  
          print(paste0(gene_output[gene, 2], ": ", nrow(df), " IVs"))
        } else {
          # gene_lookup[gene, 3] = 0
          # gene_lookup[gene, 4] = 0
          print(paste0(gene_output[gene, 2], ": Gene absent")) 
        }
      }
    }
  }
  
  write_csv(gene_output, paste0("gene_IVs_005_50_1_", colnames(ref)[path+1], ".csv"))
}