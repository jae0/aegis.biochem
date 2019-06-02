
  if (0) {
   # skip .. not ready for operational work


    if (!exists("year.assessment")) {
      year.assessment=lubridate::year(Sys.Date())
      year.assessment=lubridate::year(Sys.Date()) - 1
    }

    p = aegis::aegis_parameters( DS="biochem", yrs=1950:year.assessment )  # these are default years

   # -----------------------------
   # biochem is not yet implemented

    biochem.db( DS="biochem", p=p )

    # -----------------------------
    # stmv; vn="ca1"
    for ( vn in p$varstomodel) {
      print(vn)

      p = aegis::aegis_parameters( DS="biochem", yrs=1950:year.assessment, variables=list(Y=vn) )
      stmv( p=p, runmode=c("initialize", "globalmodel", "stage1", "stage2", "stage3", "save") )

      aegis_db( p=p, DS="predictions.redo" ) # warp predictions to other grids
      aegis_db( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
      aegis_db( p=p, DS="complete.redo" )
      aegis_db( p=p, DS="baseline.redo" )
      aegis_db_map( p=p )
      if (0) {
        global_model = stmv_db( p=p, DS="global_model")
        summary( global_model )
        plot(global_model)
      }
    }

}

}


if (0) {


      # NOTE :: integration into the stmv-approach needs to be worked out  (TODO) .. mostly done .. waiting for data to be reloaded after correctionof main table corruption


    p = list()
    p$libs = RLibrary( c( "lubridate", "lattice", "INLA", "sp" ) )
    p$libs = c( p$libs, project.library( "aegis", "stmv" )

    p$year.assessment = 2015

    p$bottles.years = c( 1955:p$year.assessment )
    p$zoop.years = c( 1990:p$year.assessment )

    recreate.aegis_5kmoffcoast=FALSE
    if (recreate.aegis_5kmoffcoast) {
      polygon_5kmoffcoast(plotdata=TRUE)
    }

    bottles.db( DS="bottles.dump.rawdata", p=p )
    bottles.db( DS="bottles.rawdata.all.redo", p=p )
    bottles.db( DS="bottles.qa.qc.redo", p=p )


    zoop.db( DS="zoop.data.dump.rawdata", p=p )
    zoop.db( DS="zoop.data.rawdata.all.redo", p=p )
    zoop.db( DS="zoop.qa.qc.redo", p=p )
    zoop.db( DS="zoop.totalSamples.redo", p=p )
    zoop.db( DS="zoop.speciesAbund.redo", p=p )


# ----


    biochem.db( DS="rawdata.datadump" )
    biochem.db( DS="flatten" )
    biochem.db( DS="scotian.shelf.redo" )



    p = list()
    p$libs = RLibrary( "lubridate", "lattice", "INLA", "sp" )
    p$libs = c( p$libs, project.library( "stmv", "aegis"  ) )


    #  biochem data system
    wd = file.path( "C:", "Gordana", "Jae", "finalData" )
    fn = "chl_buffer_filtered.csv"
    chl_fn = file.path( wd, fn)
    chl <- read.csv( chl_fn, stringsAsFactors =FALSE )
    str(chl)
    chlnames = tolower(names(chl))
    names(chl) = chlnames

    # create time stamp
    chl$ts = dmy(chl$header_start)



      # checking data
      plot( xtabs( ~ year, data=chl ))
      plot( xtabs( ~ month, data=chl ))
      plot(  header_start_lat ~  header_start_lon, data=chl, pch="." )

      i = which( chl$year %in% c(2010:2012) )
      plot(  header_start_lat ~  header_start_lon, data=chl[i,], pch="." )

      hist( chl$header_start_lat-trunc(chl$header_start_lat), "fd" )

      i = which( chl$method == 'Chl_a' )
      plot(  header_start_lat ~  header_start_lon, data=chl[i,], pch="." )
      require( lattice)
      levelplot( data_value ~ header_start_lat + chl$header_start_lon, data=chl[i,] )


      plot( data_value ~ ts, chl, pch=20, cex=.2 )

      plot( data_value ~ ts, chl[which(chl$year %in% c(1980:1984)), ], pch=20, cex=.2 )
      str(chl)



      out.dir = file.path( project.datadirectory("aegis", "biochem"), "maps", "test" )
      dir.create( out.dir, recursive=T, showWarnings=F )
      outfn = file.path( out.dir, "testmap.png" )


      p = spatial_parameters( spatial.domain="SSE" )

      datarange = seq(0,80, length.out=50)
      cols = color.code( "blue.black", datarange )

      yrs = which( chl$year %in% 2000:2005 )
      xyz = chl[, c("header_start_lon", "header_start_lat", "data_value")]
      map( xyz=xyz, depthcontours=TRUE, pts=TRUE, annot="",
           fn=outfn, loc=out.dir, at=datarange , col.regions=cols, corners=p$corners )




      #devide data into intervals
      o = cut(chl$header_start_depth, breaks =c(0,10, 30, 100,4000),include.lowest=TRUE, right=FALSE )
      m=tapply(chl$data_value,o,mean)
      s=tapply(chl$data_value,o,sd)
      n=tapply(chl$data_value,o,length)
      mm=tapply(chl$data_value,o,min)
      mx=tapply(chl$data_value,o,max)

      clim=data.frame(cbind( depth_range=names(m),mean=m,SD=s,Min=mm,Max=mx,Observations=n))
      rownames( clim) = NULL

      chl$depthcat = cut(chl$header_start_depth, breaks =c(0,10, 30, 100,4000),include.lowest=TRUE, right=FALSE )

      #tapply for multiple groups
      #m=tapply(chl$data_value, INDEX=cbind(o, chl$month) , mean)

      g = NULL
      for (m in 1:12) {
        for (d in 1:4){
          ii = which( chl$month == m & chl$depthcat == d )
          out = NULL
          out = cbind( mean=mean( chl$data_value[ii]), sd=sd(chl$data_value[ii]),
                        mm=min(chl$data_value[ii], mx=max(chl$data_value[ii]), n=length(chl$data_value[ii]))
          g = rbind( g, out)
        }
      }

      # polygon manipulations
      poly.fn = file.path(project.datadirectory("aegis", "biochem"), "datap", "BotClim_pts_v7_merge.csv"  )
      poly = read.csv( poly.fn )
      chl$polyarea = NA
      for (pl in 1:length(unique(poly$Shape_Num))) {
        pldat = poly[ which(poly$Shape_Num ==pl), ]
        test = point.in.polygon(point.x=chl$header_start_lon, point.y=chl$header_start_lat,
                                pol.x=pldat$Lon_dd, pol.y=pldat$Lat_dd )
        inside = which( test != 0 )
        chl$polyarea[inside] = pl

      }

      #other option, use tapply for multiple groups
      jj = list(o, chl$month, chl$polyarea )
      month.depth.poly =tapply(chl$data_value, INDEX=jj, mean)

      mystats = function(x) {
        m = mean(x, na.rm=TRUE)
        sd = mean(x, na.rm=TRUE)
        n =length( which(is.finite(x)))
        mx = max(x, na.rm=TRUE)
        mn = min(x, na.rm=TRUE)
        return( list( mean=m, sddev=sd, n=n, max=mx, min=mn )  )
      }


      month.depth.poly = tapply(chl$data_value, INDEX=jj, mystats)

      with(chl, tapply(chl$data_value, jj, function(x) { c(mean(x) , sd(x) )} ))




      ----
      #  space-only est of parameters
      i = which( chl$header_start_depth < 20 & chl$year == 2010 & chl$month == 6 )
      X = chl[i, c("header_start_lon", "header_start_lat", "data_value") ]
      names( X ) = c("lon","lat", "chl")
      chltest = stmv_variogram( xy=X[,c("lon", "lat")], z=X$chl , range_correlation=0.9 )



}
