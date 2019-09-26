################################################################################
# File: GSE29402_Overlap_ETL.R                                                 #
# Author: Adam Faranda                                                         #
# Created: Sept 24, 2019                                                       #
# Purpose:  load data from GSE29402 experiment data into the MySQL database    #
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
  "sample_label", "gene_id", "measurement", "meas_type"
)

###############################################################################
#                                                                             #
#                         P R E P A R E   D A T A                             #
#                                                                             #
# Import transcriptomic dataset(s) from data files or concurrent analysis.    #                                                                            
# Fix any value errors that could cause problems for MySQL                    #
###############################################################################

data_dir<-paste(wd, "data_files", sep="/")               # Path To Data File(s)

# Import Normalized Expression Values for GSE29402
setwd(data_dir)
gse<-getGEO(filename='GSE29402_series_matrix.txt.gz', destdir=data_dir)
setwd(wd)

# Import Homologs
hf<-paste(data_dir,"Mouse_Human_Homologs_MGI_090419.tsv", sep="/")
hom<-read.table(hf,  sep="\t", header=F, stringsAsFactors = F)
names(hom)[grep("V1", names(hom))]<-'Gene.Symbol'
names(hom)[grep("V5", names(hom))]<-'MGI.symbol'

dupHuman<-unique(hom[duplicated(hom$Gene.Symbol), "Gene.Symbol"])
dupMouse<-unique(hom[duplicated(hom$MGI.symbol), "MGI.symbol"])
homUnique<-hom %>% filter(!(Gene.Symbol %in% dupHuman) & !(MGI.symbol %in% dupMouse))

# Use limma to calculate differential expression statistics
print(apply(exprs(gse), 2, summary))
dm<-model.matrix(~0+source_name_ch1, pData(gse))
colnames(dm)<-gsub('source_name_ch1', '', colnames(dm))
cmat<-makeContrasts(Conjunctiva - Cornea, levels=dm)
fit<-lmFit(gse, dm)
fit<-contrasts.fit(fit, cmat)
fit<-eBayes(fit)

cols<-c(
  'ID', 'Gene.Symbol', 'ENTREZ_GENE_ID', 'AveExpr',
  'logFC', 'P.Value', 'adj.P.Val'
)
deg<-topTable(fit, n=Inf)[, cols]
print(head(deg))

# Join Group Averages on data frame with DEG statistics
COR<-pData(gse)[pData(gse)$source_name_ch1 == 'Cornea', 'geo_accession']
CNJ<-pData(gse)[pData(gse)$source_name_ch1 == 'Conjunctiva', 'geo_accession']

ex<-as.data.frame(exprs(gse))
ex$ID<-row.names(ex)
ex<-melt(ex, id.vars='ID')

gr.avg<-as.data.frame(
  inner_join(
    ex %>%
      filter(variable %in% COR) %>%
      group_by(ID) %>%
      summarize(group_1_avg=mean(value)),
    ex %>%
      filter(variable %in% CNJ) %>%
      group_by(ID) %>%
      summarize(group_2_avg=mean(value)),
    by='ID'
  )
)
deg<-as.data.frame(inner_join(deg, gr.avg,by='ID'))



# Drop Probes with ambiguos annotations (more than one gene)
deg <- deg %>% filter(!grepl('///', Gene.Symbol))

ndeg<-nrow(deg)
deg<-left_join(
  deg,
  homUnique[,c('Gene.Symbol', 'MGI.symbol')],
  by='Gene.Symbol'
)

print(
  paste(
    "Num Rows in deg before join: ", ndeg,
    ", after join:",nrow(deg),
    "With one to one murine homolog:",
    nrow(deg[!is.na(deg$MGI.symbol),])
  )
)

# Calculate minimum expression threshold per original authors
minExp<-log(0.01*(2^max(deg$AveExpr)),2)

# Format "deg" for mysql load
names(deg)[grep("^MGI.symbol$", names(deg))]<-"gene_id"
names(deg)[grep("^P.Value$", names(deg))]<-"p_value"
names(deg)[grep("^adj.P.Val$", names(deg))]<-"fdr"
deg$contrast_label<-"CORvsCNJ"



# For genes with multiple detections Unique Max Significant Fold Change
deg<-uniqueMaxLfc(deg, idc="Gene.Symbol", fdr = "fdr")

# Prepare expression values for upload
measurement<-ex %>% 
  filter(ID %in% deg$ID)

measurement<-inner_join(
  measurement, 
  deg[!is.na(deg$gene_id),c("gene_id", "ID")],
  by="ID"
)

names(measurement)[grep("^variable$", names(measurement))]<-"sample_label"
names(measurement)[grep("^value$", names(measurement))]<-"measurement"
measurement$meas_type<-"RMA Normalized Probe Intensity"


# Filter "deg" to remove human genes with no mouse homologs
pairwise<-deg[!is.na(deg$gene_id),pr_cols]


# Build Sample Table
sample <- data.frame(             # Data frame with information about samples
  label=c(CNJ, COR),
  external_accession = c(CNJ, COR),
  platform = pData(gse)[c(CNJ, COR), "platform_id"],
  proc_batch = rep(1, times=6),
  Tissue = gsub(
    "tissue: ", "",
    pData(gse)[c(CNJ, COR), "characteristics_ch1"] 
  ),
  Sample_title = pData(gse)[c(CNJ, COR), "title"],
  Minimum_Expression = rep(minExp, 6),
  stringsAsFactors = FALSE
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
  expLab = "Ramirez-Miranda",
  expTitle="Comparative gene expression in the human conjunctiva and cornea",
  expAcc = "GSE29402", 
  repo="Gene Expression Omnibus (GEO)",        
  orgn="mouse",        
  prme="RMA Normalized probe intensities, measured on an Affymetrix HG-U133_Plus 2 array (GPL570)", 
  gr_title = c(
    COR = "Corneal tissue",       
    CNJ = "Conjunctival tissue"
  ),
  gr_criteria = c(
    COR = "Tissue = cornea", 
    CNJ = "Tissue = conjunctiva"
  ),
  gr_assign = list(
    COR = COR,   
    CNJ = CNJ
  ),                           
  cn_label = c (
    CORvsCNJ = "Cornea_vs_Conjunctiva"
  ),
  cn_title = c (            # Named vector of descriptive contrast titles
    CORvsCNJ = "Cornea_vs_Conjunctiva"
  ),
  cn_sep = "vs",                       # Separator between group labels
  attr_def=c(                          # Named vector of attribute definitions
    Tissue = "the tissue providing this sample",
    Sample_title = "The descriptive title assigned to this sample by original investigator",
    Minimum_Expression = ("The minimum probe intensitty for a gene to be considered biologically significant")
  ),
  attr_units=c(
    Minimum_Expression = "Unitless"
  ),
  sampleTable = sample,
  pairwise = pairwise,
  measurement = measurement[,ms_cols]
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



