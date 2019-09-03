################################################################################
# File: Excel_Reports.R                                                        #
# Author: Adam Faranda                                                         #
# Created: June 28, 2019                                                       #
# Purpose:                                                                     #
#         Functions that write overlapping differential expression results     # 
#         to spreadsheets                                                      #
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

pairwise<-c(
  Epithelium = 'epi',
  Fibers = 'fib',
  Cornea = 'cor'
)


cm<-combn(pairwise, 2)
comparisons<-list()
for(i in 1:ncol(cm)){
  comparisons[[paste(cm[,i], collapse="_")]]<-cm[,i]
}

for(c in names(comparisons)){
  C1<-names(pairwise)[grep(comparisons[[c]][1], pairwise)]
  C2<-names(pairwise)[grep(comparisons[[c]][2], pairwise)]
  
  dg1<-ss_master %>% filter(Experiment == comparisons[[c]][1])
  dg2<-ss_master %>% filter(Experiment == comparisons[[c]][2])
  print(
    paste(
      "Comparison", c, " -- ",
      comparisons[[c]][1],"Rows:", nrow(dg1),
      comparisons[[c]][2],"Rows:", nrow(dg2)
    )
  )
  allResults<-query(dg1, dg2)
  bioResults<-bioSig(allResults)
  print(nrow(allResults))
  print(nrow(bioResults))
  
  stat.tables<-subsetTables(
    Contrast_1 = C1, Contrast_2 = C2,
    allResults, annot = an, unlog=T, descname = T
  )
  stat.inx<-tabulateOverlap(stat.tables, rename = T)
  
  bio.tables<-subsetTables(
    Contrast_1 = C1, Contrast_2 = C2,
    bioResults, annot = an, unlog=T, descname = T, stat = F
  )
  bio.inx<-tabulateOverlap(bio.tables, rename = T)
  
  print(stat.inx)
  print(bio.inx)
  
  # Set up list
  descTables = list(
    C1=list(
      Table=degSummary(dg1), corner=c(30, 2), cn=T, rn=T, tc=F
    ),
    C2=list(
      Table=degSummary(dg2), corner=c(35, 2), cn=T, rn=T, tc=F
    ),
    `Statistically Significant Intersection`=list(
      Table=stat.inx, corner=c(41, 2), cn=T, rn=T, tc=T
    ),
    `Biologically Significant Intersection`=list(
      Table=bio.inx, corner=c(45, 2), cn=T, rn=T, tc=T
    )
  )
  names(descTables)[grep("C1", names(descTables))]<-C1
  names(descTables)[grep("C2", names(descTables))]<-C2
  
  # Write Description Tables -- see script "Excel_Write_Functions.R"
  wb<-writeDescTables(
    template="Comparisons.xlsx",    # Name of template file
    name = "Data Description",      # Name of sheet with descriptive tables
    descTables = descTables
  )
  
  # Delete "Contrasts tab from directional subsets
  stat.tables["Contrasts"]<-NULL
  bio.tables["Contrasts"]<-NULL
  
  # Add statistically significant directional subsets to workbook object
  for(i in names(stat.tables)){
    writeSheet(
      wb, stat.tables[[i]], i,
      nm_cols=c("Change", "Avg")
    )
  }
  
  for(i in names(bio.tables)){
    writeSheet(
      wb, bio.tables[[i]], i,
      nm_cols=c("Change", "Avg")
    )
  }
  # Save workbooks to file  
  fn<-paste(names(pairwise)[pairwise %in% comparisons[[c]]], collapse="_")
  fn<-paste("Pax6",fn, "DEG_Comparison.xlsx", sep="_")
  #Only uncomment the following line to re-generate spreadsheets!!!!
  saveWorkbook(wb, file=fn, overwrite = T)
}
