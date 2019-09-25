################################################################################
# File: PCO_DBI_Overlap_ETL.R                                                  #
# Author: Adam Faranda                                                         #
# Created: June 27, 2019                                                       #
# Purpose: Template script for loading experiment data into the MySQL database #
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
dna<-list()
for(f in list.files(data_dir, pattern="WT0vs")){
  c<-gsub(".xlsx","", f)
  file<-paste(data_dir,f, sep="/")
  dna[[c]]<-read.xlsx(file)[,-c(1,2)]
  names(dna[[c]])[grep('^gene$', names(dna[[c]]))]<-'gene_id'
  names(dna[[c]])[grep('q_value', names(dna[[c]]))]<-'fdr'
  names(dna[[c]])[grep('sample_1', names(dna[[c]]))]<-'Group_1'
  names(dna[[c]])[grep('sample_2', names(dna[[c]]))]<-'Group_2'
  names(dna[[c]])[grep('value_1', names(dna[[c]]))]<-'group_1_avg'
  names(dna[[c]])[grep('value_2', names(dna[[c]]))]<-'group_2_avg'
  
  dna[[c]]$Group_1<-gsub("W_h0", "WT0", dna[[c]]$Group_1)
  dna[[c]]$Group_2<-gsub(
    "W_h6", "WT6",  
    gsub("W_h24", "WT24", dna[[c]]$Group_2)
  )
  
  dna[[c]]$contrast_label<-paste(
    dna[[c]]$Group_1, 
    dna[[c]]$Group_2,
    sep='vs'
  )
  
  dna[[c]]$Platform<-"RNASeq"
  dna[[c]]$Interval<-gsub(
    "WT6", "6_Hour",
    gsub(
      "WT24", "24_Hour",
      gsub(
        "WT48", "48_Hour", dna[[c]]$Group_2
      )
    )
  )
  
  dna[[c]]$Lab<-"DNA"
  dna[[c]]$Experiment<-"PCO"
  dna[[c]]$logFC<-log2(dna[[c]]$group_2_avg/dna[[c]]$group_1_avg)
  dna[[c]]<-dna[[c]] %>% filter(status=="OK")
}

# Build Sample Table
sample <- data.frame(             # Data frame with information about samples
  label=paste("S", 1:3, sep=""),
  external_accession = rep(NA, times=3),
  platform = rep("RNA-Seq", 3),
  proc_batch = rep(1, times=3),
  Genotype = rep("WT", times=3),  
  Age = rep(5, times=3),
  Interval = rep(c(0, 6, 24), each=1),
  stringsAsFactors = FALSE
)

pairwise<-data.frame(stringsAsFactors = F)
for(c in names(dna)){
  pairwise<-bind_rows(
    pairwise,
    dna[[c]][,pr_cols]
  )
}

measurement = data.frame(
  sample_label = character(0),
  gene_id = character(0),
  measurement = numeric(0),
  meas_type = character(0)
)

# Verify that attribute gene_id in tables for measurements and pairwise results
# contains values <= 50 characters in length
pairwise$gene_id<-strtrim(pairwise$gene_id, 50)

# Replace any "Inf" or "-Inf" logFC values with an appropriate, easily
# recognized numerical value 
print(
  paste(
    "Maximum Absolute Foldchange in this dataset is:", 
    max(abs(pairwise$logFC))
  )
)
infval<-1000
print(paste("Replacing Inf with:", infval))
pairwise$logFC[pairwise$logFC == Inf]<-infval
pairwise$logFC[pairwise$logFC == -Inf]<- 0-infval

# Verify that gene identifiers are unique within each contrast or sample
for(c in unique(pairwise$contrast_label)){
  gid<-pairwise[pairwise$contrast_label == c,"gene_id"]
  print(
    paste(
      "Rows in contrast",c,":", nrow(pairwise[pairwise$contrast_label == c,]),
      "Unique gene_id:", length(unique(gid))
    )
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
  con=cnx,                
  expLab = "DNA Link",
  expTitle="mouse lens epithelial cells after cateract surgery ",
  expAcc = "None", 
  repo="Not Public",        
  orgn="mouse",        
  prme="TaKaRa Pico Paired End Sequencing, Tuxedo Suite Analysis", 
  gr_title = c(
    WT0 = "Wildtype, zero hours PCS",       
    WT6 = "Wildtype, 6 hours PCS",     
    WT24 = "Wildtype, 24 hours PCS"         
  ),
  gr_criteria = c(
    WT0 = "Genotype = WT, Interval = 0", 
    WT6 = "Genotype = WT, Interval = 6",
    WT24 = "Genotype = WT, Interval = 24"
  ),
  gr_assign = list(
    WT0 = c("S1"),   
    WT6 = c("S2"),   
    WT24 = c("S3")   
  ),                           
  cn_label = c (
    WT0vsWT6 = "Wildtype_0_vs_6_Hours", 
    WT0vsWT24 = "Wildtype_0_vs_24_Hours"
  ),
  cn_title = c (            # Named vector of descriptive contrast titles
    WT0vsWT6 = "LEC from Wildtype mice: 0 hours PCS vs 6 Hours PCS",
    WT0vsWT24 = "LEC from Wildtype mice: 0 hours PCS vs 24 Hours PCS"
  ),
  cn_sep = "vs",                     # Separator between group labels
  attr_def=c(                          # Named vector of attribute definitions
    Genotype = "Sample genotype WT = Wildtype",
    Age = "The age of mice at the time of the experiment",
    Interval = ("Time elapsed after removal of lens fibers prior to capsule harvest")
  ),
  attr_units=c(
    Age = "Months",
    Interval = "Hours"
  ),
  sampleTable = sample,
  pairwise = pairwise,
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



