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

# Example Call to load an ensembl GTF File
# gtfp <- "/home/adam/Documents/LTS_Data/Mus_musculus.GRCm38.96.chr.gtf"
# etl_load_gtf_table(
#   gtfPath = gtfp,
#   con=cnx,
#   gtf_table_name = "ENSMUS_96_GTF"
# )
# etl_build_gene_table_ens(con=cnx)

# Load a Marker Table from MGI -- Note
mgi_file<-"/home/adam/Desktop/MGI_Gene_Model_Coord.rpt"
etl_load_gene_table_mgi(
  con = cnx,                                  # the database connection to use
  mgi_file = mgi_file,                        # Name of the data file to load
  mgi_table_name = "MGI",                     # Name of the table to load into
  idCol = 'MGI_accession_id',                 # Column with unique identifier
  syCol = 'marker_symbol',                    # Column with unique gene symbols
  dsCol = 'marker_name',                      # Column with descriptive name
  colType = 'VARCHAR(500)'                    # Column type for fields in MySQL
)
# Note on MGI Columns idCol, syCol, dsCol
# In the original file, the names of these columns have a whitespace which
# is replaced with an underscore in the above function call.  In the body
# of the function, there is a line that changes the column headers after
# loading data into a data frame.  The values of idCol, syCol and dsCol

