###############################################################################
# File: Gene_Table_ETL.R                                                      #
# Purpose: Populate the database table GENE with gene information extracted   #
#          from a GTF file                                                    #
# Created: August 10, 2019                                                    #
# Author: Adam Pater-Faranda                                                  #
###############################################################################

options(stringsAsFactors = F)
source('DEG_Overlap_ETL_Functions.R')
wd<-getwd()
loaddir<-"/var/lib/mysql-files/"  # Directory where "txt" data must be saved
library(dplyr)
library(data.table)
library(DBI)
library(RMySQL)
library(rtracklayer)
library(biomaRt)
cnx<-dbConnect(
  MySQL(), 
  dbname="DEG_Overlap_Analysis_System", 
  password="adnaraf7"
)

# Push GTF File to MySQL database
gtfp <- "/home/adam/Documents/LTS_Data/Mus_musculus.GRCm38.96.chr.gtf"
etl_load_data(
  df=etl_load_gtf_table(con=cnx, gtfPath = gtfp),
  con=cnx,
  table_name = "ENSMUS_96_GTF",
  prefix=""
)

# Retrieve Gene Symbols that correspond to 
etl_load_data(
  df = etl_build_gene_table_ens(con=cnx),
  con = cnx,
  table_name = "GENE"
)

