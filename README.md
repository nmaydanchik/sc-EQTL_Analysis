# Multivariable Mendelian randomization analysis mapping risk genes with cell-type-specific effects on Alzheimer’s Disease
This repository contains coding samples from work I did for a biostatistics laboratory at UChicago in an Educational Assignment where I used multivariable Mendelian Randomization methods to analyze single-cell expression Quantitative Trait Loci data to determine **the causal effects of cell-type-specific gene expression on Alzheimer's Disease**. 

Project Dates: Jan-Feb 2024 (Junior year of High School)

This repository highlights my ability to:

* Process, reconcile, combine, and modify datasets with millions of values in R
* Conduct Mendelian randomization biostatistical analysis
* Learn new advanced topics (genomics, biostatistics) and apply it in analysis
* Generate plots (ggplot2)
* Create presentations

I presented my findings to a UChicago lab group of ~15 Ph.D. students and postdocs as well as two professors 2/19/24. The presentation is available as FebPresentation.pdf. This presentation consists of two projects; the work shared in this repository relates to the second project.

## Project and File Description

The objective of this project was to map genes with cell-type-specific effects on Alzheimer's Disease (focusing on cells in brain cortex tissues). This was part of a greater research effort of a Ph.D. student that culminated in a published paper in The American Journal of Human Genetics.

First, I did a focused analysis of the CAV1 gene in endothelial cells. The scripts involved in this analysis are CAV1_Analysis.R, as well as supporting scripts RunMethods.R and process_scEQTL_data.R.

These were the instructions given to me by the Ph.D. student:

Input:
* Alzheimer's disease GWAS summary statistics
* sc-eQTL summary statistics, compare with bulk tissue eQTL summary statistics  

Analysis:
* Filter sc-eQTL summary statistics and only keep CAV1 gene.
* Do LD clumping (250kb, r^2<0.3, or other criteria), P<0.005.
* Perform MR methods.  

Expected output:
* Cell-type-specific effects of CAV1 on AD.

I found that the LD clumping parameters and the IV-gene expression effect p-value restriction were too stringent and resulted in not enough instrumental variables (IV) for analysis (1, to be precise). I softened the restrictions to 250kb, r^2<0.5, p<0.05, and indeed found a consistent negative effect of expression of the CAV1 gene in endothelial cells on Alzheimer's Disease. 

Afterwards, I did a multivariable Mendelian randomization analysis where I analyzed nearly 2,000 gene/cell-type combinations, split into seven biological pathways for Alzheimer's Disease, looking for genes with cell-type-specific effects on Alzheimer's disease. The script involved in this analysis is mv_scEQTL_analysis.R. The findings of this analysis are in the presentation.
