################################################################################
# File: DEG_Plotting_Functions.R                                               #
# Author: Adam Faranda                                                         #
# Created: Sept 24, 2019                                                       #
# Purpose:  Plotting functions for analyzing overlap                           #
#                                                                              #
################################################################################
# Setup Workspace
rm(list=ls()) # Clear any pre-existing environment data
options(echo=F)
library('openxlsx')
library('dplyr')
library('reshape2')
library('GEOquery')
library('limma')
wd<-getwd()
source('Overlap_Comparison_Functions.R')
source('DEG_Overlap_ETL_Functions.R')

gene_id_target = "symbol"

# Standardized Column Headers for pairwise result table
pr_cols<-c(
  "contrast_label", "gene_id", "group_1_avg", "group_2_avg","logFC", 
  "p_value", "fdr"
)

# Standardized Column Headers for measurement table
ms_cols<-c(
  "sample_label", "gene_id", "measurement_value", "meas_type"
)

