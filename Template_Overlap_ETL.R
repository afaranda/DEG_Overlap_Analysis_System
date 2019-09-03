################################################################################
# File: Template_Overlap_ETL.R					                                       #
# Author: Adam Faranda			                                     				       #
# Created: June 27, 2019						                                           #
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

gene_id_target = "external_id"

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

# Example Sample Data table 
sample <- data.frame(             # Data frame with information about samples
  label=c("S1", "S2", "S3", "S4", "S5", "S6"),
  external_accession = c(
    'GSM123', 'GSM134', 'GSM135',
    'GSM223', 'GSM234', 'GSM235'
  ),
  platform = rep("RNA-Seq",6),
  proc_batch = c(1,1,1,1,1,1),
  Genotype = rep(c("WT", "FN", "P6"), each=2),  
  Age = c(1,10,100, 1, 10, 100),
  Gender = c('M', 'M', 'F', 'M', 'M', 'F'),
  stringsAsFactors = FALSE
)

# Exaple Pairwise Result table
pairwise <- data.frame(
  contrast_label = c("10yr_WTvsFN", "10yr_FNvsP6"),
  gene_id = c("ENSMUSG00000000003", "ENSMUSG00000000003"),
  group_1_avg = c(8, 4),
  group_2_avg = c(4, 8),
  logFC = c(-1, 1),
  p_value = c(0.05, 0.0001),
  fdr = c(0.5, 0.01),
  stringsAsFactors = F
)

# Example Sample Measurement table
measurement = data.frame(
  sample_label = c("S1", "S2", "S3", "S4", "S5", "S6"),
  gene_id = rep("ENSMUSG00000000003", times = 6),
  measurement = c(6, 10, 2, 6, 6, 10),
  meas_type = rep("count", times = 6)
)
# Verify that attribute gene_id in tables for measurements and pairwise results
# contains values <= 50 characters in length


# Replace any "Inf" or "-Inf" logFC values with an appropriate, easily
# recognized numerical value 

# Verify that gene identifiers are unique within each contrast or sample

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
  con=cnx,                 # the database connection to use
  expLab = "The Lab",  # The lab where this experiment was
  expTitle="title",    # The descriptive title of this experiment
  expAcc = "SRA12345", # The external accession number for this experiment
  repo="SRA",          # The source of the raw data
  orgn="mouse",        # The organism this data set represents
  prme="the only way", # The method used to process sampes in this dataset
  gr_title = c(
    Group_1 = "The control Group",       # Named vector of descriptive group
    Group_2 = "The Treatment Group",     # names, (Names ARE group
    Group_3 = "The Other Group"          # labels)
  ),
  gr_criteria = c(
    Group_1 = "Genotype = WT, Age < 10", # Named vector of grouping criteria
    Group_2 = "Genotype = FN, Age < 10", # names, (Names MUST match group
    Group_3 = "Genotype = P6, Age < 10"  # labels)
  ),
  gr_assign = list(
    Group_1 = c("S1", "S2"),   # List of Group assignments, each
    Group_2 = c("S3", "S4"),   # Element is named by group label
    Group_3 = c("S5", "S6")    # and contains a vector of
  ),                           # sample labels
  cn_label = c (
    Group_1_vs_Group_2 = "10yr_WTvsFN",  # Named vector of Labels for contrasts
    Group_2_vs_Group_3 = "10yr_FNvsP6",  # Names must match group labels,
    Group_1_vs_Group_3 = "10yr_WTvsP6"   # separated by value of 'cn_sep'
  ),
  cn_title = c (            # Named vector of descriptive contrast titles
    Group_1_vs_Group_2 = "Aged mice Wildtype vs Pax6 Null",
    Group_2_vs_Group_3 = "Aged mice Fibr Null vs Pax6 Null",
    Group_1_vs_Group_3 = "Aged mice Wildtype vs Fibr Null"
  ),
  cn_sep = "_vs_",                     # Separator between group labels
  attr_def=c(                          # Named vector of attribute definitions
    Genotype = "Sample genotype WT = Wildtype, P6 = Pax6 Null, FN = FN Null",
    Age = ("Age of the mouse at the time of sample collection")
  ),
  attr_units=c(                        # Named vector of attribute units
    Age = "Years"
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



