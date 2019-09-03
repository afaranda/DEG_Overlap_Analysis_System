################################################################################
# File: Generic_Excel_Report.R                                                 #
# Author: Adam Faranda                                                         #
# Created: June 28, 2019                                                       #
# Purpose:                                                                     #
#         Generic Reporting script used to create overlap spreadsheets for     #
#         one or more pairs of pairwise contrasts.  This script should be      #
#         updated to reflect the desired query                                 #
################################################################################

# Setup Workspace   
options(echo=T)
library('openxlsx')
library('dplyr')
library('org.Mm.eg.db')
source("Pax6Epi_Multi_Compare_Overlap.R")
source("Excel_Write_Functions.R")
load('master_tables.Rdata')


# Build annotation table  -- just use the AnnotationDbi Package directly.
# This annotation table does not need to change, as long as it is keyed
# to the genes in the query set
an<-data.frame(
  AnnotationDbi::select(
    org.Mm.eg.db, 
    keys = unique(ag_master$MGI.symbol),
    columns = c("SYMBOL", "GENENAME"),
    keytype = "SYMBOL"
  ) %>% 
    group_by(SYMBOL) %>%
    filter(row_number() == 1),
  stringsAsFactors = F
)
names(an)[grep("SYMBOL", names(an))]<-"MGI.symbol"
names(an)[grep("GENENAME", names(an))]<-"description"

# Defines a set of pairwise contrasts: Element names
# are used in the output, element values correspond
# to specific values in ss_master / ag_master
pairwise<-c(
  DNA_Link = 'DNA',
  DBI = 'DBI'
)

# Constructs a list of pairs of comparisons to process
# Each list element
cm<-combn(pairwise, 2)
comparisons<-list()
for(i in 1:ncol(cm)){
  comparisons[[paste(cm[,i], collapse="_")]]<-cm[,i]
}

for(c in names(comparisons)){
  for(c in names(comparisons)){
    
    # Extract Names from comparison list object
    C1<-names(pairwise)[grep(comparisons[[c]][1], pairwise)]
    C2<-names(pairwise)[grep(comparisons[[c]][2], pairwise)]
    
    # Prepare DEG sets for the desired comparison.  Each deg set
    # MUST be keyed on a unique Gene Symbol / ID etc for the
    # subsequent join.
    dg1<-ss_master %>% filter(
      Lab == comparisons[[c]][1] & Group_2 == "WT24"
    )
    dg2<-ss_master %>% filter(
      Lab == comparisons[[c]][2] & Group_2 == "WT24"
    )
    
    # Validate Key in Joined tables ()
    print(
      paste(
        "Comparison", c, " -- ",
        comparisons[[c]][1],"Rows:", nrow(dg1),
        "Unique Genes:", length(unique(dg1$MGI.symbol)),
        comparisons[[c]][2],"Rows:", nrow(dg2),
        "Unique Genes:", length(unique(dg2$MGI.symbol))
      )
    )
    
    # Builds a spreadsheet for the comparison 
    createBioSigOverlapSpreadSheet(
      C1 = C1, C2 = C2, template ="Comparisons.xlsx",
      dg1 = dg1, dg2 = dg2, pref = "PCO_24", 
      dg1.ds ="Samples Analyzed at", dg2.ds = "Samples Analyzed at"
    )
  }
}
