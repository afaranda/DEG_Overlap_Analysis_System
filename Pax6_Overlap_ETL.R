################################################################################
# File: Pax6_Overlap_ETL.R		    			                                       #
# Author: Adam Faranda			                                     				       #
# Created: June 27, 2019						                                           #
# Purpose:  #
#                                                                              #
################################################################################
# Setup Workspace
rm(list=ls()) # Clear any pre-existing environment data
options(echo=F)
library('openxlsx')
library('dplyr')
library('reshape2')

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

###############################################################################
#                                                                             #
#                         P R E P A R E   D A T A                             #
#                                                                             #
# Import transcriptomic dataset(s) from data files or concurrent analysis.    #                                                                            
# Fix any value errors that could cause problems for MySQL                    #
###############################################################################

data_dir<-paste(wd, "data_files", sep="/")               # Path To Data File(s)

#Import Corneal WT vs Pax6 Het pairwise contrast 
corDEG<-read.xlsx(paste(data_dir, 'Pax6_Cornea_DEG.xlsx', sep="/"))
corDEG<-corDEG[,setdiff(names(corDEG), c("gene_id", "test_id"))]

# Fix Column Headers for Cornea Pax6DE
names(corDEG)[grep('^gene$', names(corDEG))]<-'gene_id'
names(corDEG)[grep('q_value', names(corDEG))]<-'fdr'
names(corDEG)[grep('sample_2', names(corDEG))]<-'Group_1'
names(corDEG)[grep('sample_1', names(corDEG))]<-'Group_2'
names(corDEG)[grep('value_2', names(corDEG))]<-'group_1_avg'
names(corDEG)[grep('value_1', names(corDEG))]<-'group_2_avg'
corDEG$contrast_label<-"cor_WTvsP6"
corDEG$logFC<-log2(corDEG$group_1_avg/corDEG$group_2_avg)
corDEG<-corDEG %>% filter(status=="OK")

# Import Lens Epithelial Cell WT vs Pax6 Het pairwise contrast 
epiDEG<-read.xlsx(paste(data_dir, 'Cuffdiff1_max.xlsx', sep="/"))
epiDEG<-epiDEG[,setdiff(names(epiDEG), c("gene_id", "test_id"))]

# Fix Column Headers for Lens Epithelium Pax6 DEG
names(epiDEG)[grep('^gene$', names(epiDEG))]<-'gene_id'
names(epiDEG)[grep('q_value', names(epiDEG))]<-'fdr'
names(epiDEG)[grep('sample_1', names(epiDEG))]<-'Group_1'
names(epiDEG)[grep('sample_2', names(epiDEG))]<-'Group_2'
names(epiDEG)[grep('value_1', names(epiDEG))]<-'group_1_avg'
names(epiDEG)[grep('value_2', names(epiDEG))]<-'group_2_avg'
epiDEG$contrast_label<-"epi_WTvsP6"
epiDEG$logFC<-log2(epiDEG$group_1_avg/epiDEG$group_2_avg)
epiDEG<-epiDEG %>% filter(status=="OK")

# Import Lens Fiber CEll WT vs Pax6 Het pairwise contrast
fibDEG<-read.xlsx(paste(data_dir,'Cuffdiff2_max.xlsx', sep="/"))
fibDEG<-fibDEG[,setdiff(names(fibDEG), c("gene_id", "test_id"))]

# Fix Column Headers for Lens fibers Pax6 DEG
names(fibDEG)[grep('^gene$', names(fibDEG))]<-'gene_id'
names(fibDEG)[grep('q_value', names(fibDEG))]<-'fdr'
names(fibDEG)[grep('sample_1', names(fibDEG))]<-'Group_1'
names(fibDEG)[grep('sample_2', names(fibDEG))]<-'Group_2'
names(fibDEG)[grep('value_1', names(fibDEG))]<-'group_1_avg'
names(fibDEG)[grep('value_2', names(fibDEG))]<-'group_2_avg'
fibDEG$contrast_label<-"fib_WTvsP6"
fibDEG$logFC<-log2(fibDEG$group_1_avg/fibDEG$group_2_avg)
fibDEG<-fibDEG %>% filter(status=="OK")

# Combine Pairwise Results into a single table
master<-bind_rows(
  epiDEG[,pr_cols],
  fibDEG[,pr_cols],
  corDEG[,pr_cols]
)

# There is one record in each of the DNA Link data sets where the 
# gene symbol is a long string of comma separated "Igkv-XXX" values 
# (~1551 characters).  The following line fixes this error by trimming
# values in "ag_master$gene_id" to a maximum of 50 characters
master$gene_id<-strtrim(master$gene_id, 50)

# MySQL requires actual values for numeric columns
# Replace logFC values of "Inf" or "-Inf" with 1000 or -1000
# The actual maximum observed (Non-infite) fold change is 11.6
master$logFC[master$logFC == Inf]<-1000
master$logFC[master$logFC == -Inf]<- -1000


# Enter information about samples
sample <- data.frame(            
  label=c("S1", "S2", "S3", "S4", "S5", "S6"),
  external_accession = c(NA, NA, NA, NA, NA, NA),
  platform = rep("RNA-Seq", 6),
  proc_batch = c(1,1,2,2,2,2),
  Genotype = rep(c("WT", "P6"), times=3),  
  Age = c(5,5,5,5,5,5),
  Tissue = c('cor', 'cor', 'epi', 'epi', 'fib', 'fib'),
  stringsAsFactors = FALSE
)

# Sample Measurement table -- No measurement data is being loaded with this set
measurement = data.frame(
  sample_label = character(),
  gene_id = character(),
  measurement = character(),
  meas_type = character()
)
# Verify that gene identifiers are unique within each contrast or sample
x<-list(
  master %>% filter(contrast_label == "epi_WTvsP6"),
  master %>% filter(contrast_label == "fib_WTvsP6"),
  master %>% filter(contrast_label == "cor_WTvsP6")
)
for (i in 1:length(x)){
  print(
    paste("Rows in ",i, ":", nrow(x[[i]]), 
          "Unique gene_id: ", length(unique(x[[i]]$gene_id)))
  )
}

###############################################################################
#                                                                             #
#                   L O A D   S T A G I N G   T A B L E S                     #
#                                                                             #
#    Edit the arguments to the function 'etl_load_experiment_tables'          #
#    such that they reflect the actual experiment                             #                                                                            
#                                                                             #
###############################################################################

# delete staging tables from any previous loads
dbSendQuery(cnx, "CALL drop_staging_tables();")

# Load data and meta data into staging tables
etl_load_experiment_tables(
  con = cnx,
  expLab="DNA Link",
  expTitle ="Pax6 dependent genes",
  expAcc = "None",
  repo = "Not Public",
  orgn = "Mouse",
  prme = "TaKaRa Pico Paired End Sequencing, Tuxedo Suite Analysis",
  gr_title = c(
    cor_WT = "Pax6 Wildtype Homozygous Cornea",        
    cor_P6 = "Pax6 Null/Wildtype Heterozygous Cornea",
    epi_WT = "Pax6 Wildtype Homozygous Lens Epithelium",        
    epi_P6= "Pax6 Null/Wildtype Heterozygous Lens Epithelium",
    fib_WT = "Pax6 Wildtype Homozygous Lens Fibers",        
    fib_P6= "Pax6 Null/Wildtype Heterozygous Lens Fibers"
  ),
  gr_criteria = c(
    cor_WT= "Genotype = 'WT' AND Tissue = 'cor'",
    cor_P6 = "Genotype = 'P6' AND Tissue = 'cor'",
    epi_WT = "Genotype = 'WT' AND Tissue = 'epi'",
    epi_P6 = "Genotype = 'P6' AND Tissue = 'epi'",
    fib_WT = "Genotype = 'WT' AND Tissue = 'fib'",
    fib_P6  = "Genotype = 'P6' AND Tissue = 'fib'"
  ),
  gr_assign = list( 
    cor_WT = c("S1"), cor_P6 = c("S2"),   
    epi_WT = c("S3"), epi_P6 = c("S4"),
    fib_WT = c("S5"), fib_P6 = c("S6")
  ),                      
  cn_label = c (
    cor_WT_vs_cor_P6 = "cor_WTvsP6",
    epi_WT_vs_epi_P6 = "epi_WTvsP6",
    fib_WT_vs_fib_P6= "fib_WTvsP6"
  ),
  cn_title = c (            # Named vector of descriptive contrast titles
    cor_WT_vs_cor_P6 = "Wildtype vs Pax6 Heterozygous Cornea",
    epi_WT_vs_epi_P6 = "Wildtype vs Pax6 Heterozygous Lens Epithelium",
    fib_WT_vs_fib_P6 = "Wildtype vs Pax6 Heterozygous Lens Fibers"
  ),
  cn_sep = "_vs_",                     # Separator between group labels
  attr_def=c(                          # Named vector of attribute definitions
    Genotype = "Sample genotype WT = Wildtype, P6 = Pax6 Heterozygous",
    Age = "Age of the mouse at the time of sample collection",
    Tissue = "The Tissue from which this sample was collected"
  ),
  attr_units=c(                        # Named vector of attribute units
    Age = "Months"
  ),
  sampleTable = sample,
  pairwise = master,
  measurement = measurement
)


###############################################################################
#                                                                             #
#                   M I G R A T E   A N D    R E P O R T                      #
#                                                                             #
#    migrate data from staging tables into main system, report the number of  #
#    rows loaded for each table.                                              #                                                                            
#                                                                             #
###############################################################################

# Migrate data from staging tables into main database tables
query<-paste(
  "CALL load_pairwise_tables('", gene_id_target, "', @lr);", sep=""
)
dbSendQuery(conn=cnx, query)
res<-dbSendQuery(conn=cnx, "SELECT @lr;")
print(dbFetch(res))

# report number of genes skipped when loading pairwise tables
print("Number of genes dropped from pairwise results")
print(
  dbFetch(
    dbSendQuery(
      conn=cnx,
      "SELECT COUNT(*) FROM vw_pairwise_nogene;"
    )
  )
)

# report number of genes skipped when loading sample measurement tables
print("Number of genes dropped from sample measurements")
print(
  dbFetch(
    dbSendQuery(
      conn=cnx,
      "SELECT COUNT(*) FROM vw_pairwise_nogene;"
    )
  )
)

# Summary Count of Records loaded in this transaction
print("Number of records loaded into each table for this transaction")
print(
  dbFetch(
    dbSendQuery(
      conn=cnx,
      "CALL report_counts()"
    )
  )
)



