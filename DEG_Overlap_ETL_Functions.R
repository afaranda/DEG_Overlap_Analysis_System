###############################################################################
# File: DEG_Overlap_ETL_Script.R                                              #
# Purpose: Load data from a differential expression experiment into tables    #
#          in the DEG Overlap database. This script needs to be run as root   #
# Created: August 10, 2019                                                    #
# Author: Adam Pater-Faranda                                                  #
###############################################################################

# Setup workspace
options(stringsAsFactors = F)
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


###############################################################################
# Function: etl_load_data                                                     # 
# Given a data frame, where column names correspond to                        #
# attribute names in the target database table, export the data frame as tab  #
# separated text, build and run a "LOAD DATA INFILE" query, then delete the   #
# data file after confirming the load                                         #
###############################################################################
etl_load_data<-function(
  df,                   # data frame to load into the database
  con,                  # conection object (points to target database)
  prefix="",            # prefix for temporary filename
  table_name="SAMPLE",   # name of target table in database
  id_null=F             # Whether or not to append "SET id = NULL" to query
){
  fn<-paste(loaddir, prefix, "_", table_name, "_etl.txt", sep="")
  write.table(
    df, file=fn, row.names = F, col.names = F, quote = T, sep=","
  )  
  
  
  query<-paste("'", fn, "'", sep="")
  query<-paste("LOAD DATA INFILE", query)
  query<-paste(query, "INTO TABLE", table_name)
  query<-paste(query, "FIELDS TERMINATED BY ','")
  query<-paste(query, "OPTIONALLY ENCLOSED BY '\"'")
  colnames<-paste("(", paste( names(df), collapse=", "), ")" )
  if(id_null){
    query<-paste(query, colnames, "SET id = NULL;")
  } else {
    query <-paste(query, " ",colnames, ";", sep="")
  }
  print(query)
  
  result<-NULL
  result<-dbSendQuery(con, query)
  if(!is.null(result)){
    if(file.exists(fn)){file.remove(fn)}
  }
}
###############################################################################
# Function: etl_load_gtf_table()                                              #
# load a gtf file into its own table on the database                          #
#                                                                             #
###############################################################################
gtfp <- "/home/adam/Documents/LTS_Data/Mus_musculus.GRCm38.96.chr.gtf"
etl_load_gtf_table<-function(
  con,                                # MySQL database connection   
  gtfPath = gtfp,                     # Path to gtf file
  gtf_table_name = "ENSMUS_96_GTF",   # Name of target table in MySQL database
  colType = "VARCHAR(200)"            # Column type for fields in MySQL table
){
  gtf<-readGFF(
    gtfPath
  )
  gtfCols<-paste(
    paste(names(gtf), colType),  collapse =", "
  )
  create_statement<-paste(
    "(id INT NOT NULL AUTO_INCREMENT,", gtfCols, ", PRIMARY KEY (id));"
  )
  
  query<-paste(
    "CREATE TABLE IF NOT EXISTS", gtf_table_name, 
    create_statement
  )
  dbSendQuery(con, query)
  
  etl_load_data(
    df = gtf,
    con = con, 
    table_name = gtf_table_name 
  )
}

###############################################################################
# Function: etl_build_gene_table_ens                                          #
#           Transfer gene information from an ensmbl gtf file into the table  #
#           'GENES'. Add gene specifiec information (descriptions)            #
#           retrieved from ensembl via biomaRt.  This function returns a      #
#           data frame suitable for loading to the database via the function  #
#           'etl_load_data'                                                   #
###############################################################################

etl_build_gene_table_ens<-function(
  con = cnx,
  gtfTable = 'ENSMUS_96_GTF',
  idCol = 'gene_id', 
  nmCol = 'gene_name',
  tpCol = "`type`",
  ds = "mmusculus_gene_ensembl",
  attr = c(
    "ensembl_gene_id",   
    "mgi_symbol", 
    "mgi_description"
  ), 
  filt = c("ensembl_gene_id")
){
  min_attr<-c(
    "ensembl_gene_id",   
    "mgi_symbol", 
    "mgi_description"
  )
  if(any(!(min_attr %in% attr))) {
    print("attribute set must contain the following: ")
    print("       'ensembl_gene_id', 'mgi_symbol, and 'mgi_description")
    return(NULL)
  } else {
    # Query the gtf table to get id's and symbols that match 1 to 1
    query<-paste("SELECT GROUP_CONCAT(", idCol,") AS ", idCol,",", sep="")
    query<-paste(query, "gene_name, COUNT(*) AS nrow")
    query<-paste(query,"FROM", gtfTable, "WHERE",tpCol,"='gene'")
    query<-paste(query,"GROUP BY",nmCol,"HAVING nrow < 2 LIMIT 100000;")
    rs<-dbSendQuery(con, query)
    gn<-dbFetch(rs, n=-1)
    
    print(paste("Number of 1 to 1 genes:", nrow(gn)))
    # Fetch Descriptions from Ensembl via biomaRt for genes in "gn"
    an<-getBM(
      attributes = attr,
      filters = filt,
      values = unique(gn[,idCol]),
      mart = useMart(
        "ensembl",
        dataset = "mmusculus_gene_ensembl"
      )
    )
    print(paste("Number of genes retrieved from biomaRt:", nrow(an)))
    
    # Merge Query genes, fix column headers
    df<-merge(
      x=gn,
      y=an,
      by.x=ifelse(is.numeric(idCol), names(gn)[idCol], idCol),
      by.y='ensembl_gene_id'
    )
    print(paste("Rows in df after merge:", nrow(df)))
    
    # Filter for genes where the MGI symbol matches the gtf gene_name
    df<-df[df$gene_name == df$mgi_symbol, ]
    print(paste("Rows in df after filtering symbols:", nrow(df)))
    
    # Update Column Names, add "external source"
    ex_id<-ifelse(is.numeric(idCol), names(gn)[idCol], idCol)
    names(df)[grep(nmCol, names(df))]<-"symbol"
    names(df)[grep(ex_id, names(df))]<-"external_id"
    names(df)[grep("mgi_description", names(df))]<-"description"
    df$external_id_source<-"ensembl"
    
    # Setup 'df' for loading
    cols<-c('symbol', 'external_id', 'external_id_source', 'description')
    df<-df[,cols]
    etl_load_data(
      df = df,
      con = con,
      table_name = "GENE"
    )
  }
}
###############################################################################
# Function: etl_load_gene_table_mgi                                           #
#   Transfer gene information from an mgi marker file into the table 'GENES'. #
#                                                                             #
###############################################################################
etl_load_gene_table_mgi <- function(
  con,                                        # the database connection to use
  mgi_file = "MGI_Gene_Model_Coord.rpt",      # Name of the data file to load
  mgi_table_name = "mgi",                     # Name of the table to load into
  idCol = 'MGI_accession_id',                 # Column with unique identifier
  syCol = 'marker_symbol',                    # Column with unique gene symbols
  dsCol = 'marker_name',                      # Column with descriptive name
  colType = 'VARCHAR(500)'                    # Column type for fields in MySQL
){
  df<-read.table(
    file=mgi_file, header=T, sep="\t",
    quote="", comment.char = ""
  )
  print(nrow(df))
  # Make column headers sql compliant
  names(df)<-gsub("\\.", "_", names(df))
  
  # Verify that at least one identifying column is
  if(
      !(
        length(unique(df[,idCol])) == nrow(df) |
        length(unique(df[,syCol])) == nrow(df)
      ) 
  ){
    print("Table must provide either Unique Identifiers or Unique Symbols")
    return(NULL)
  }
  if(any(names(df) == "id")){
    print("The header 'id' is reserved; please edit this header in the file")
    return(NULL)
  }
  
  
  # Create table on the database that matches the MGI File
  tblCols<-paste(
    paste(names(df), colType),  collapse =", "
  )
  create_statement<-paste(
    "(id INT NOT NULL AUTO_INCREMENT,", tblCols, ", PRIMARY KEY (id));"
  )
  query<-paste(
    "CREATE TABLE IF NOT EXISTS", mgi_table_name, 
    create_statement
  )
  dbSendQuery(con, query)
  
  # Push the MGI File into the database
  etl_load_data(
    df=df,
    con=con,
    table_name = mgi_table_name
  )
  
  # Update Column Names for downstream filtering
  names(df)[grep(paste("^",idCol,"$", sep=""), names(df))]<-"idCol"
  names(df)[grep(paste("^",syCol,"$", sep=""), names(df))]<-"syCol"
  names(df)[grep(paste("^",dsCol,"$", sep=""), names(df))]<-"dsCol"
  
  
  # Remove genes without 1 to 1 correspondence between MGI symbol and ID
  df<-df  %>% 
    dplyr::select(idCol, syCol, dsCol)
  
  dg <- as.data.frame(
    inner_join(df, df, by = "syCol") %>% 
      group_by(idCol.x) %>%
      filter(n() < 2) %>%
      dplyr::select(idCol=idCol.x, syCol, dsCol=dsCol.x)
  )
  
  # Prepare Dataframe for loading into the database
  names(dg)[grep("idCol", names(dg))]<-"external_id"
  names(dg)[grep("syCol", names(dg))]<-"symbol"
  names(dg)[grep("dsCol", names(dg))]<-"description"
  dg$external_id_source<-"MGI"
  
  print(names(dg))
  etl_load_data(
    df=dg,
    con=con, 
    table_name = "GENE"
  )
}

###############################################################################
# Function: etl_load_pairwise_tables                                          #
# Given an custom gene list and some additional information, divide           #
# data into separate data frames and load them into their corresponding       #
# database tables                                                             #
###############################################################################
etl_load_experiment_tables<-function(
  con,                 # the database connection to use
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
  ),                                     # sample labels
  cn_label = c (
    Group_1_vs_Group_2 = "10yr_WTvsFN",  # Named vector of Labels for contrasts
    Group_2_vs_Group_3 = "10yr_FNvsP6",  # Names must match group labels, 
    Group_1_vs_Group_3 = "10yr_WTvsP6"   # separated by "_vs_"
  ),
  cn_title = c (            # Named vector of descriptive contrast titles
    Group_1_vs_Group_2 = "Aged mice Wildtype vs Pax6 Null",  
    Group_2_vs_Group_3 = "Aged mice Fibr Null vs Pax6 Null", 
    Group_1_vs_Group_3 = "Aged mice Wildtype vs Fibr Null"  
  ),
  cn_sep = "_vs_",          # Separator between group labels
  sampleTable = data.frame(             # Data frame with information about samples
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
  ),
  attr_def=c( # Named vector of attribute definitions
    Genotype = "Sample genotype WT = Wildtype, P6 = Pax6 Null, FN = FN Null",
    Age = ("Age of the mouse at the time of sample collection")
  ),
  attr_units=c( # Named vector of attribute units
    Age = "Years"
  ),
  pairwise = data.frame(
    contrast_label = c("10yr_WTvsFN", "10yr_FNvsP6"),
    gene_id = c("ENSMUSG00000000003", "ENSMUSG00000000003"),
    group_1_avg = c(8, 4),
    group_2_avg = c(4, 8),
    logFC = c(-1, 1),
    p_value = c(0.05, 0.0001),
    fdr = c(0.5, 0.01),
    stringsAsFactors = F
  ),
  measurement = data.frame(
    sample_label = c("S1", "S2", "S3", "S4", "S5", "S6"),
    gene_id = rep("ENSMUSG00000000003", times = 6),
    measurement = c(6, 10, 2, 6, 6, 10),
    meas_type = rep("count", times = 6)
  )
){
  # Validate Data Being Pushed to the database
  pc_dup<-as.data.frame(
    pairwise %>% 
      group_by(contrast_label, gene_id) %>%
      summarize(Duplicates=n()) %>%
      filter(Duplicates > 1)
  )
  if(nrow(pc_dup) > 0){
    print("Duplicate genes detetcted in pairwise contrast table")
    print("First 10 Duplicates:")
    print(pc_dup[,1:10])
    return(pc_dup)
  }
  
  ms_dup<-as.data.frame(
    measurement %>% 
      group_by(sample_label, gene_id) %>%
      summarize(Duplicates=n()) %>%
      filter(Duplicates > 1)
  )
  if(nrow(ms_dup) > 0){
    print("Duplicate genes detetcted in sample measurement table")
    print("First 10 Duplicates:")
    print(ms_dup[,1:10])
    return(ms_dup)
  }
  # Create Staging Tables
  dbSendQuery(con, "CALL create_staging_tables();")
  #Load Table stg_EXPERIMENT
  etl_load_data(
    df = data.frame(
      lab = expLab,
      title = expTitle,
      repository = repo,
      external_accession = expAcc,
      organism = orgn,
      proc_method = prme,
      stringsAsFactors = F
    ), con = con, table_name = "stg_EXPERIMENT", id_null = T
  )
  # Load Table SAMPLE
  if(any(duplicated(sampleTable$label))){
    print("All Samples must be labeled uniquely")
    return (NULL)
  } else {
    st_cols<-intersect(
      names(sampleTable),
      dbListFields(con, "stg_SAMPLE")
    )
    at_cols<-setdiff(
      names(sampleTable),
      dbListFields(con, "stg_SAMPLE")
    )
    etl_load_data(
      df = sampleTable[,st_cols], con = con, table_name = "stg_SAMPLE"
    )
    
    # Load Tables SAMPLE_NUM_ATR and SAMPLE_CAT_ATR
    at_cols<-sapply(sampleTable[,at_cols], class)
    
    num<-data.table::melt(sampleTable, id.vars='label', variables.factor=F) %>%
      filter(variable %in% names(at_cols[at_cols == "numeric"]))
    num$variable<-as.character(num$variable)
    num$attr_def<-sapply(num$variable, function(x) attr_def[x])
    num$attr_units<-sapply(num$variable, function(x) attr_units[x])
    
    names(num)[grep("label", names(num))]<-"sample_label"
    names(num)[grep("variable", names(num))]<-"attr_name"
    names(num)[grep("value", names(num))]<-"attr_value"
    
    etl_load_data(
      df=num,
      con = con,
      table_name = "stg_SAMPLE_NUM_ATR"
    )
    
    chr<-data.table::melt( sampleTable, id.vars='label', variables.factor=F) %>%
      filter(variable %in% names(at_cols[at_cols != "numeric"]))
    chr$variable<-as.character(chr$variable)
    chr$attr_def<-sapply(chr$variable, function(x) attr_def[x])
    
    names(chr)[grep("label", names(chr))]<-"sample_label"
    names(chr)[grep("variable", names(chr))]<-"attr_name"
    names(chr)[grep("value", names(chr))]<-"attr_value"
    
    etl_load_data(
      df = chr,
      con = con,
      table_name = "stg_SAMPLE_CAT_ATR"
    )
  }
  # Load Table stg_GROUP
  if(
    any(names(gr_title) != names(gr_criteria)) |
    any(names(gr_title) != names(gr_assign)) 
  ){
    print("Group label mismatch")
    print("Group labels must match, and must be in the same order")
    return(NULL)
  } else {
    groups<-names(gr_title)
    grp<-data.frame(
      label = groups,
      title = sapply(groups, function(x) gr_title[x]),
      criteria = sapply(groups, function(x) gr_criteria[x]),
      stringsAsFactors = F
    )
    etl_load_data(
      df = grp, 
      con = con,
      table_name = "stg_GROUP" 
    )
    
    # Load Table stg_GROUP_ASSIGN
    gr_as<-melt(gr_assign)
    gr_as$value<-as.character(gr_as$value)
    names(gr_as)[grep("value", names(gr_as))]<-"sample_label"
    names(gr_as)[grep("L1", names(gr_as))]<-"group_label"
    
    etl_load_data(
      df =  gr_as,
      con = con,
      table_name = "stg_GROUP_ASSIGN"
    )
  }
  # Load Table stg_CONTRAST
  if(any(names(cn_label) != names(cn_title))){
    print("Contrast label mismatch")
    print("Contrast labels must match, and must be in the same order")
  } else {
    cn_as<-data.frame(
      label=cn_label,
      title=cn_title[names(cn_label)],
      group_1_label=sapply(
        strsplit(names(cn_label), cn_sep),
        function(x) x[1]
      ),
      group_2_label=sapply(
        strsplit(names(cn_label), cn_sep),
        function(x) x[2]
      ),
      row.names = 1:length(cn_label),
      stringsAsFactors = F
    )
    etl_load_data(
      df=cn_as,
      con=con,
      table_name = "stg_CONTRAST"
    )
  }
    
  # Load Table stg_PAIRWISE_RESULT
  etl_load_data(
    df=pairwise,
    con = con, 
    table_name = 'stg_PAIRWISE_RESULT'
  )
  
  # Load Table stg_SAMPLE_MEASUREMENT
  etl_load_data(
    df=measurement,
    con = con, 
    table_name = 'stg_SAMPLE_MEASUREMENT'
  )
  
}

###############################################################################
# Function: etl_load_dgelist                                                  #
# Given a "dge_list" object and some additional information, divide data into #
# separate data frames and load them into their corresponding database tables #
###############################################################################

###############################################################################
# Function: etl_load_eset                                                     #
# Given an "expressionSet" object and some additional information, divide     #
# data into separate data frames and load them into their corresponding       #
# database tables                                                             #
###############################################################################



# pair_ctr = data.frame(
#   contrast = rep(c("10yr_WTvsFN", "10yr_FNvsP6"), each = 10),
#   genes_id = rep(paste("ENMUSG00000123", 0:9, 0:9, sep=""), times=2),
#   group_1_avg = rep(8, times=20),
#   group_2_avg = rep(c(4, 16), times=10),
#   logFC = rep(c(-1, 1), times=10),
#   p_value = rep(c(0.05, 0.0001),times=10),
#   fdr = rep(c(0.5, 0.01),times=10),
#   stringsAsFactors = F
# )



