################################################################################
# File: PCO_DNA_Overlap_ETL.R                                                  #
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

dbi<-list()
for(f in list.files(data_dir, pattern="_expressedTags-all.txt")){
  c<-gsub("_expressedTags-all.txt","", f)
  c<-gsub("_vs_", "vs", c)
  c<-gsub("WT_0_Hour", "WT0", c)
  c<-gsub("WT_24_Hour", "WT24", c)
  c<-gsub("WT_48_Hour", "WT48", c)
  
  file<-paste(data_dir,f, sep="/")
  dbi[[c]]<-read.table(
    file, sep="\t",
    header= T, quote="", stringsAsFactors = F
  )
  samples<-names(dbi[[c]])[grep("_RPKM", names(dbi[[c]]))]
  print(paste("RPKM Filtering Cols:", paste(samples, collapse=" "), sep=" "))
  dbi[[c]]$group_1_avg<-apply(dbi[[c]][,samples[1:3]], 1, mean)
  dbi[[c]]$group_2_avg<-apply(dbi[[c]][,samples[4:6]], 1, mean)
  
  dbi[[c]]$Group_1<-unlist(strsplit(c,"vs"))[1]
  dbi[[c]]$Group_2<-unlist(strsplit(c,"vs"))[2]
  dbi[[c]]$contrast_label<-paste(
    dbi[[c]]$Group_1, 
    dbi[[c]]$Group_2,
    sep='vs'
  )
  
  dbi[[c]]$Platform<-"RNASeq"
  dbi[[c]]$Interval<-gsub(
    "WT6", "6_Hour",
    gsub(
      "WT24", "24_Hour",
      gsub(
        "WT48", "48_Hour", dbi[[c]]$Group_2
      )
    )
  )
  
  names(dbi[[c]])[grep('gene_id$', names(dbi[[c]]))]<-'ensembl'
  names(dbi[[c]])[grep('GeneID$', names(dbi[[c]]))]<-'gene_id'
  names(dbi[[c]])[grep('PValue$', names(dbi[[c]]))]<-'p_value'
  names(dbi[[c]])[grep("FDR", names(dbi[[c]]))]<-'fdr'
  
  names(dbi[[c]])<-gsub("_GeneCount.counts$", "", names(dbi[[c]]))
  dbi[[c]]$Lab<-"DBI"
  dbi[[c]]$Experiment<-"PCO"
  dbi[[c]]<-uniqueMaxLfc(dbi[[c]], idc='gene_id', fdr='fdr')
}

# Build Sample Table
sample <- data.frame(             # Data frame with information about samples
  label=paste("S", 1:9, sep=""),
  external_accession = rep(NA, times=9),
  platform = rep("RNA-Seq", 9),
  proc_batch = rep(1, times=9),
  Genotype = rep("WT", times=9),  
  Age = rep(5, times=9),
  Interval = rep(c(0, 24, 48), each=3),
  stringsAsFactors = FALSE
)

pairwise<-data.frame()
for(c in names(dbi)){
  pairwise<-bind_rows(
    pairwise,
    dbi[[c]][,pr_cols]
  )
}

measurement<-data.frame()
for(c in names(dbi)){
  cols<-c("gene_id", intersect(names(dbi[[c]]), sample$label))
  measurement<-rbind(
    measurement, melt(dbi[[c]][,cols]))
}

# Since Samples 1 - 3 are control in both cases, they show up twice.  Verify that they are 
# Identical in both cases
check<-measurement %>% 
  arrange(gene_id, variable) %>% 
  group_by(gene_id, variable) %>%
  summarise(Val1 = nth(value, 1), Val2 = nth(value, 2), diff =  nth(value, 1) - nth(value, 2)) %>%
  filter(!is.na(Val2))
print(paste("All duplicate control values for S1 - S3 are identical:", !any(check$diff != 0)))


# Format table for upload
names(measurement)[grep('variable', names(measurement))]<-"sample_label"
names(measurement)[grep('value', names(measurement))]<-"measurement"
measurement$meas_type <-"Count"
measurement <- as.data.frame(
  measurement %>%
    group_by(gene_id, sample_label) %>%
    filter(duplicated(gene_id) == F)
)
# Verify that attribute gene_id in tables for measurements and pairwise results
# contains values <= 50 characters in length
pairwise$gene_id<-strtrim(pairwise$gene_id, 50)
measurement$gene_id<-strtrim(measurement$gene_id, 50)

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

for(c in unique(measurement$sample_label)){
  gid<-measurement[measurement$sample_label == c,"gene_id"]
  print(
    paste(
      "Rows for this sample",c,":", nrow(measurement[measurement$sample_label == c,]),
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
  con=cnx,                 # the database connection to use
  expLab = "DBI",
  expTitle="mouse lens epithelial cells after cateract surgery ",
  expAcc = "None", 
  repo="Not Public",        
  orgn="mouse",        
  prme="TaKaRa Pico Paired End Sequencing, Tuxedo Suite Analysis", 
  gr_title = c(
    WT0 = "Wildtype, zero hours PCS",       
    WT24 = "Wildtype, 24 hours PCS",     
    WT48 = "Wildtype, 48 hours PCS"         
  ),
  gr_criteria = c(
    WT0 = "Genotype = WT, Interval = 0", 
    WT24 = "Genotype = WT, Interval = 24",
    WT48 = "Genotype = WT, Interval = 48"
  ),
  gr_assign = list(
    WT0 = c("S1", "S2", "S3"),   
    WT24 = c("S4", "S5", "S6"),   
    WT48 = c("S7", "S8", "S9")   
  ),                           
  cn_label = c (
    WT0vsWT24 = "Wildtype_0_vs_24_Hours", 
    WT0vsWT48 = "Wildtype_0_vs_48_Hours"
  ),
  cn_title = c (            # Named vector of descriptive contrast titles
    WT0vsWT24 = "LEC from Wildtype mice: 0 hours PCS vs 24 Hours PCS",
    WT0vsWT48 = "LEC from Wildtype mice: 0 hours PCS vs 48 Hours PCS"
  ),
  cn_sep = "vs",                       # Separator between group labels
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



