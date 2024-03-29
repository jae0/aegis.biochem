

biochem_parameters = function( p=NULL, project_name=NULL, project_class="core", ... ) {


  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args


  # ---------------------
  # create/update library list
  p$libs = c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "GADMTools" ) )
  p$libs = c( p$libs, project.library ( "aegis", "aegis.biochem" ) )

  p$project_class = project_class

  p$project_name = ifelse ( !is.null(project_name), project_name, "biochem" )

  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project_name )
  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=F, recursive=T )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=F, recursive=T )

  if (!exists("spatial_domain", p) ) p$spatial_domain = "SSE"
  if (!exists("spatial_domain_subareas", p)) p$spatial_domain_subareas = c( "snowcrab", "SSE.mpa" )
  p = spatial_parameters( p=p)

  # define focal years
  if (!exists( "yrs", p)) p$yrs = 1950:lubridate::year(lubridate::now())  # default
  p = temporal_parameters(p=p, aegis_dimensionality="space-year" )


  if (project_class=="core") {
    return(p)
  }


  if (project_class %in% c( "stmv") ) {
    p$libs = c( p$libs, project.library ( "stmv" ) )
    p$DATA = 'biochem.db( p=p, DS="stmv_inputs" )'
    p$varstomodel = c( )
    if (!exists("stmv_variables", p)) p$stmv_variables = list()
    if (!exists("LOCS", p$stmv_variables)) p$stmv_variables$LOCS=c("plon", "plat")
    if (!exists("TIME", p$stmv_variables)) p$stmv_variables$TIME="tiyr"
    p = aegis_parameters(p=p, DS="stmv") # generics
    return(p)
  }



  if (project_class %in% c( "carstm" ) ) {
    p$libs = c( p$libs, project.library ( "carstm" ) )
    return(p)
  }



  if (project_class %in% c( "hybrid", "default") ) {
    p$libs = c( p$libs, project.library ( "carstm" ) )
    return(p)
  }


}
