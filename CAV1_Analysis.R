library(tidyverse)
library(data.table)
library(ieugwasr)
library(MendelianRandomization)
library(mr.raps)
library(MRMix)
library(GRAPPLE)
library(patchwork)
setwd("C://Users/Nate/Documents/RLibrary")
source("RunMethods.R")
source("process_sc_eQTL_Data.R")


# Process Endothial Cell Data (Bryois) ------------------------------------

s <- proc.time()[3]
X <- fread("sc_eQTL_AllCells/EndothialCells/Endothelial.cells.1.gz", select = c(1,2,4,5))
# SNP_Alleles <- fread("snp_pos.txt.gz", select = c("SNP", "effect_allele", "other_allele")) 

for(chromnum in 2:22) {
  filename = paste0("sc_eQTL_AllCells/EndothialCells/Endothelial.cells.",chromnum,".gz")
  X <- rbind(X, fread(filename, select = c(1,2,4,5)))
  print(paste("Chromosome", chromnum, "done"))
}
f <- proc.time()[3]
print(f-s)

s <- proc.time()[3]
filenames = list.files("sc_eQTL_AllCells/EndothialCells", pattern = "*.gz", full.names = TRUE)
X <- lapply(filenames, fread, select = c(1,2,4,5)) %>% bind_rows()
f <- proc.time()[3]
print(f-s)

colnames(X) = c("Gene_id", "rsid", "pval", "beta")
X = X[grep('CAV1', X$Gene_id)]

X <- mutate(X, se = abs(X$beta/qnorm(X$pval/2)))

Y <- fread("ieu_AD.txt.gz") %>% 
  select(!"chr")
colnames(Y) = c("non_effect_allele.ad", "effect_allele.ad", "beta.ad", "se.ad", "pval.ad", "rsid")

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
df <- select(df, !effect_allele:effect_allele.ad)

filtered_df <- filter(df, pval < 0.05) %>% 
  ld_clump(clump_kb = 250, clump_r2 = 0.5)

mr.obj = mr_input(bx = filtered_df$beta, bxse = filtered_df$se, by = filtered_df$beta.ad, byse = filtered_df$se.ad)
est <- RunMethods(filtered_df$rsid, mr.obj, allMethods = FALSE)
res <- mr_ivw(mr.obj)

write_csv(as.data.frame(est), "mr_res.csv")

geneinfo <- read_tsv("gene_info_GTEx_v8.txt.zip")

# Find number of IVs for cell type ------------------------------------------------------------

res = matrix(nrow = 14, ncol = 2)
colnames(res) = c("cell_type", "numIV")

# Load Alleles corresponding to eQTL data
SNP_Alleles <- fread("snp_pos.txt.gz", select = c("SNP", "effect_allele", "other_allele")) 

# Load Alzheimer's Disease data
Y <- fread("ieu_AD.txt.gz") %>% 
  select(!"chr")
colnames(Y) = c("non_effect_allele.ad", "effect_allele.ad", "beta.ad", "se.ad", "pval.ad", "rsid")

for (filenum in 1:14) {
  print(paste0("file: ", filenum))
  df <- process_sc_eQTL_Data(filenum, study = "Jerber")
  if (!is.null(df)) {
    res[filenum, 1] = df$cellClusterName[1]
    res[filenum, 2] = nrow(df)    
  }
}

# Run MR on specific cell type ---------------------------------------------

df <- processBryoisData(6) # Get data
mr.obj = mr_input(bx = df$beta, bxse = df$se, by = df$beta.ad, byse = df$se.ad)
est <- RunMethods(df$rsid, mr.obj)

# Bulk Tissue eQTL analysis ----------------------------------------------

# Load eQTL Data
X <- fread("Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz", select = c(1,2,7,8,9))
X = X[grep('ENSG00000105974', X$gene_id)]
lookup = fread("GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz", select = c(1,4,5,7))

# Merge eQTL data with reference info to get rsid and alleles
X <- inner_join(X, lookup, by = join_by(variant_id))
rm(lookup); gc();
X <- select(X, c(3:8))
colnames(X) = c("pval", "beta", "se", "effect_allele", "non_effect_allele", "rsid")

# Load Alzheimer's Disease data
Y <- fread("ieu_AD.txt.gz") %>% 
  select(!"chr")
colnames(Y) = c("non_effect_allele.ad", "effect_allele.ad", "beta.ad", "se.ad", "pval.ad", "rsid")

# Merge with AD data
df <- inner_join(X, Y, by = join_by(rsid))

# Match alleles between the eQTL and GWAS data
for (i in 1:nrow(df)) {
  if (df[i,]$effect_allele == df[i,]$non_effect_allele.ad && df[i,]$non_effect_allele == df[i,]$effect_allele.ad) {
    df[i,]$beta = -1 * df[i,]$beta # If they are switched, invert beta of eQTL
  }
  else if (df[i,]$effect_allele != df[i,]$effect_allele.ad || df[i,]$non_effect_allele != df[i,]$non_effect_allele.ad)  {
    print(paste0("Allele mismatch: row ", i)) # If they are not switched and not the same, output a warning
  }
}
df <- select(df, !c(effect_allele, non_effect_allele, effect_allele.ad, non_effect_allele.ad))

filtered_df <- df %>% 
  ld_clump(clump_kb = 250, clump_r2 = 0.5)

mr.obj = mr_input(bx = filtered_df$beta, bxse = filtered_df$se, by = filtered_df$beta.ad, byse = filtered_df$se.ad)
est <- RunMethods(filtered_df$rsid, mr.obj, allMethods = FALSE)
res <- mr_ivw(mr.obj)


# Plots

p1 <- df %>%
  ggplot(mapping = aes(x = beta, y = beta.ad)) +
  geom_point(alpha = 0.8, color = "#00BFC4") +
  geom_linerange(aes(xmin = beta - se, xmax = beta + se), alpha = 0.3) +
  geom_linerange(aes(ymin = beta.ad - se.ad, ymax = beta.ad + se.ad), alpha = 0.3) +
  labs(
    x = "Effect of SNP on CAV1 Expression",
    y = "Effect on SNP on Alzheimer's Disease",
    subtitle = "Slope of Regression: -0.2864",
    title = "Without LD Clumping"
  ) +
  coord_cartesian(ylim = c(-.25, .01), xlim = c(.35, 1.05)) +
  geom_smooth(method = "lm")
p2 <- filtered_df %>%
  ggplot(mapping = aes(x = beta, y = beta.ad)) +
  geom_point(alpha = 0.8, color = "#00BFC4") +
  geom_linerange(aes(xmin = beta - se, xmax = beta + se), alpha = 0.3) +
  geom_linerange(aes(ymin = beta.ad - se.ad, ymax = beta.ad + se.ad), alpha = 0.3) +
  labs(
    x = "Effect of SNP on CAV1 Expression",
    y = "Effect on SNP on Alzheimer's Disease",
    subtitle = "Slope of Regression: -0.2183",
    title = "With LD Clumping"
  ) +
  coord_cartesian(ylim = c(-.25, .01), xlim = c(.35, 1.05)) +
  geom_smooth(method = "lm")
p <- p1 + p2 + plot_annotation(
  title = "Effect sizes of SNPs on CAV1 Expression and Alzheimer's Disease",
  subtitle = "Bulk Tissue"
)
p

p1 <- filtered_df %>%
  ggplot(mapping = aes(x = beta, y = beta.ad)) +
  geom_point(alpha = 0.8, color = "#00BFC4") +
  geom_linerange(aes(xmin = beta - se, xmax = beta + se), alpha = 0.3) +
  geom_linerange(aes(ymin = beta.ad - se.ad, ymax = beta.ad + se.ad), alpha = 0.3) +
  labs(
    x = "Effect of SNP on CAV1 Expression",
    y = "Effect on SNP on Alzheimer's Disease",
    subtitle = "Slope of Regression: -0.0237"
  ) +
  coord_cartesian(ylim = c(-.08, .08)) +
  geom_smooth(method = "lm")
p2 <- filtered_df2 %>%
  ggplot(mapping = aes(x = beta, y = beta.ad)) +
  geom_point(alpha = 0.8, color = "#00BFC4") +
  geom_linerange(aes(xmin = beta - se, xmax = beta + se), alpha = 0.3) +
  geom_linerange(aes(ymin = beta.ad - se.ad, ymax = beta.ad + se.ad), alpha = 0.3) +
  labs(
    x = "Effect of SNP on CAV1 Expression",
    y = "Effect on SNP on Alzheimer's Disease",
    subtitle = "Slope of Regression: -0.0535",
    title = "Outlier Removed"
  ) +
  coord_cartesian(ylim = c(-.08, .08)) +
  geom_smooth(method = "lm")
p <- p1 + p2 + plot_annotation(
  title = "Effect sizes of SNPs on CAV1 Expression and Alzheimer's Disease",
  subtitle = "Endothial Cell Tissue. p < .05, Clump window = 250kb, R-squared < 0.05"
)
p
