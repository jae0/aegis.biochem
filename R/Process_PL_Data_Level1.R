Process_PL_Data_Level1 <- function(sql_file, data_source, username, password, output_file) {

  # Function to extract data from database, perform basic quality control and save the data in .Rdata and
  # tsv (tab seperated value) formats.
  #
  # Input: sql_file : sql filename - string
  #        data_source : data source name for database - string
  #        username : username - string
  #        password: user password - string
  #        output_file : output filename - string
  #
  # Output: None
  #
  # Last update: 20141015
  # Benoit.Casault@dfo-mpo.gc.ca

  # source custom functions
  # source("Extract_BioChem_Data.R")

  # get the data
  data <- Extract_BioChem_Data(sql_file, data_source, username, password)

  # rename variables
  sql_str <- data[[1]]
  df <- data[[2]]
  rm(data)

  # rename df variables to lower case
  col_names <- tolower(names(df))
  names(df) <- col_names

  # find data type (i.e. whether sample_id's, counts or weight data)
  samples <- FALSE
  weights <- FALSE
  counts <- FALSE
  if ("wet_weight" %in% col_names | "dry_weight" %in% col_names) {
    weights <- TRUE
  } else if ("counts" %in% col_names) {
    counts <- TRUE
  } else {
    samples <- TRUE
  }

  # quality control checks
  # checks on sample_id's
  if (samples) {
    # check depth values
    index <- df$start_depth==df$end_depth
    if (any(index)) {
      index <- which(index)
      sample_str <- paste(df$custom_sample_id[index], collapse=',')
      warning("Line 38: cases where start_depth=end_depth were found\n",
              "  Sample_ID: ", sample_str,"\n",
              "  These samples were removed from further analysis")
      # remove samples from dataframe
      df <- df[-index,]
    }
    # check volume values
    index <- df$volume==0
    if (any(index)) {
      index <- which(index)
      sample_str <- paste(df$custom_sample_id[index], collapse=',')
      warning("Line 49: cases where volume=0 were found\n",
              "  Sample_ID: ", sample_str, "\n",
              "  These samples were removed from further analysis")
      # remove samples from dataframe
      df <- df[-index,]
    }
    # checks on data
  } else if (weights | counts) {
    # check split fraction
    index <- df$split_fraction==0 | is.na(df$split_fraction)
    if (any(index)) {
      index <- which(index)
      sample_str <- paste(df$custom_sample_id[index], collapse=',')
      warning("Line 61: cases where split_fraction=0 or NA were found\n",
              "  Sample_ID: ", sample_str,"\n",
              "  These samples were removed from further analysis")
      # remove samples from dataframe
      df <- df[-index,]
    }
    # check sieve sizes  Need some work because of NA is sieve entries
    #       index <- df$min_sieve==df$max_sieve
    #       if (any(index)) {
    #         index <- which(index)
    #         sample_str <- paste(df$custom_sample_id[index], collapse=',')
    #         warning("Line 61: cases where min_sieve==max_sieve were found\n",
    #                 "  Sample_ID: ", sample_str,"\n",
    #                 "  These samples were removed from further analysis")
    #         # remove samples from dataframe
    #         df <- df[-index,]
    #       }
    # There should be another check to detect multiple values of split fraction for given parameter
  }

  # save to RData file
  save(file=output_file, list=c("df","sql_str"))

  # print to tsv (tab seperated values) file
  output_file <- paste(substr(output_file, 1, nchar(output_file)-5), "tsv", sep="")
  write.table(df, file=output_file, row.names=FALSE, na="", sep="\t")

  return()
}
