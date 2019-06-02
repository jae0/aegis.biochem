Extract_BioChem_Data <- function(sql_file, data_source, username, password) {

  # Function to extract data from a database (e.g. Biochem).  The function reads a sql query, opens a
  # connection to a database, extracts the data and returns the extracted data and the associated
  # sql code.
  #
  # Input: sql_file : sql filename - string
  #        data_source : data source name for database - string
  #        username : username - string
  #        password: user password - string
  #
  # Output: data : list of two elements
  #         data[[1]] : sql code for the query - string
  #         data[[2]] : data extracted from the database - dataframe
  #
  # Last update: 20141015
  # Benoit.Casault@dfo-mpo.gc.ca

  # source custom functions
  # source("Read_SQL.R")

  # declare empty list to store outputs
  data <- list()

  # read sql file
  data[[1]] <- Read_SQL(sql_file);

  # open connection to BioChem
  conn <- ROracle::dbConnect( DBI::dbDriver("Oracle"), dbname=data_source, username=username, password=password, believeNRows=F)

  # get data from BioChem
  data[[2]] <-  ROracle::dbGetQuery(conn, data[[1]], nullstring="")

  # close BioChem connection
  ROracle::dbDisconnect(conn)

  # clear variables
  rm(conn)

  return(data)
}
