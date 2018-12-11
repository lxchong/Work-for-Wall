####################################################################################################
# Author: Luke Chong
# Date: June 2 2018
#
# Adapted from Ivan Marin-Franch's original vflayoutmw_singleField.r code to generate voronoi tessellation
# continous probability plot of a single field.
####################################################################################################
require(deldir)

####################################################
# extracts x, y coordinates from list of tiles and adds
# dist to the distance of each point from the origin
####################################################

edgePntFn <- function(tile, dist = 5) {
  x <- tile$pt['x']
  y <- tile$pt['y']
  len <- sqrt(x^2 + y^2) + dist
  ang <- atan2(y, x)
  return(c(x = len * cos(ang), y = len * sin(ang)))
}

#####################################################################################################
# Need to define new vfplot function
#
#####################################################################################################

vfplot <- function (vf, plotType, xmin = NULL, xmax = NULL, ymin = NULL, 
          ymax = NULL, notSeenAsBlack = TRUE, newWindow = FALSE, txtfont = "sans", 
          pointsize = 10, width = 6, showaxis = FALSE, colaxis = "white") {
  
  if (nrow(vf) > 1) {
    stop("Error! vf cannot have more than 1 rows")
  }
  bsxy <- NULL
  bsxy$x <- c(15, 15)
  bsxy$y <- c(3, -3)
  bsxy <- as.data.frame(bsxy)
  locini <- vfsettings$locini
  evaltxt <- paste("vfsettings$", vf$tpattern, "$locnum", sep = "")
  loc_num <- eval(parse(text = evaltxt))
  evaltxt <- paste(vf$tperimetry, "locmap$", vf$tpattern, sep = "")
  patternMap <- eval(parse(text = evaltxt))
  patternMap <- patternMap[, c("xod", "yod")]
  xrange <- max(patternMap$xod) - min(patternMap$xod)
  yrange <- max(patternMap$yod) - min(patternMap$yod)
  if (is.null(xmin)) 
    xmin <- min(patternMap$xod) - 0.025 * xrange
  if (is.null(xmax)) 
    xmax <- max(patternMap$xod) + 0.025 * xrange
  if (is.null(ymin)) 
    ymin <- min(patternMap$yod) - 0.025 * yrange
  if (is.null(ymax)) 
    ymax <- max(patternMap$yod) + 0.025 * yrange
  if (plotType == "vf") {
    dev <- vf
  }
  if (plotType == "td") {
    dev <- tdval(vf)
    devp <- tdpmap(dev)
  }
  if (plotType == "pd") {
    dev <- pdval(tdval(vf))
    devp <- pdpmap(dev)
  }
  if (plotType == "pdghr") {
    dev <- pdvalghr(tdval(vf))
    devp <- pdpmapghr(dev)
  }
  if (plotType == "vf") {
    plotColor <- vfgrayscale(dev[, locini:(locini + loc_num - 
                                             1)], age = vf$sage, pattern = vf$tpattern, algorithm = vf$talgorithm)
    cloneDev <- as.character(round(dev[, locini:(locini + 
                                                   loc_num - 1)]))
    cloneDev[which(dev[, locini:(locini + loc_num - 1)] < 
                     0)] = "<0"
  }
  else {
    plotColor <- vfcolormap(as.numeric(devp[, locini:(locini + 
                                                        loc_num - 1)]))
    cloneDev <- as.numeric(dev[, locini:(locini + loc_num - 
                                           1)])
    cloneDev <- as.character(round(cloneDev))
    if (notSeenAsBlack) {
      idxblack <- which(vf[locini:(locini + loc_num - 1)] <= 
                          0)
      if (length(idxblack) > 0) 
        plotColor[idxblack, ] <- 0
    }
  }
  plotColor[is.na(plotColor)] <- 0
  for (i in 1:nrow(bsxy)) {
    if (xmin < bsxy$x[i] & xmax > bsxy$x[i] & ymin < bsxy$y[i] & 
        ymax > bsxy$y[i]) {
      idx <- which(patternMap$xod == bsxy$x[i] & patternMap$yod == 
                     bsxy$y[i])
      if (length(idx) > 0) {
        if (plotType == "vf") {
          plotColor[idx, ] <- c(0, 0, 0)
        }
        else {
          plotColor[idx, ] <- c(0.5, 0.5, 0.5)
        }
      }
      else {
        patternMap <- rbind(patternMap, c(bsxy$x[i], 
                                          bsxy$y[i]))
        cloneDev <- c(cloneDev, NA)
        if (plotType == "vf") {
          plotColor <- rbind(plotColor, c(0, 0, 0))
        }
        else {
          plotColor <- rbind(plotColor, c(0.5, 0.5, 0.5))
        }
      }
    }
  }
  if (vf$seye == "OS") {
    xmin2 <- xmin
    xmin <- -xmax
    xmax <- -xmin2
    patternMap$xod <- -patternMap$xod
  }
  vftess <- vftessellation(patternMap, dist = 3)
  vftiles <- tile.list(vftess[[1]])
  vfhull <- vftess[[2]]
  height <- width * (ymax - ymin)/(xmax - xmin)
  if (newWindow) {
    device <- options("device")
    if (.Platform$OS.type == "unix") {
      if (Sys.info()["sysname"] == "Darwin") {
        options(device = "quartz")
        dev.new(width = width, height = height, dpi = 85)
      }
      else {
        options(device = "x11")
        dev.new(width = width, height = height)
      }
    }
    else {
      options(device = "windows")
      dev.new(width = width, height = height, rescale = "fixed")
    }
    options(device = device)
  }
  vfplotloc(cloneDev, patternMap = patternMap, vftiles = vftiles, 
            vfhull = vfhull, loccol = plotColor, xmin = xmin, xmax = xmax, 
            ymin = ymin, ymax = ymax, txtfont = txtfont, pointsize = pointsize, 
            showaxis = showaxis, colaxis = colaxis)
 
  #create white convex hull for printout LXC 07/06/2018
  if (any(vf$tpattern == c("pPeri_v2","pPeriv_v2","pPerivi_v2","pPeriv"))) {
    if (vf$seye == "OD") {
      vtess <- deldir(saplocmap$pPC26$xod,saplocmap$pPC26$yod)
    } else {
      vtess <- deldir(saplocmap$pPC26$xod * -1,saplocmap$pPC26$yod)
    }
    vtessTL <- tile.list(vtess)
    outerEdges <- as.data.frame(t(sapply(vtessTL, edgePntFn, dist = 2)))
    outerEdges <- outerEdges[chull(outerEdges),]
    names(outerEdges) <- c('x', 'y')
    polygon(outerEdges,col = "white",border = "white")
  }
  
  #Add horizontal and vertical midlines
  axis( 1, pos = 0, labels = FALSE, lwd.ticks = 0, at = c( xmin - 10, xmax + 10 ), col = "black" )
  axis( 2, pos = 0, labels = FALSE, lwd.ticks = 0, at = c( ymin - 10, ymax + 10 ), col = "black" )
}

########################################################################################################################
# Function to plot voronoi tessellation printout
#
#
########################################################################################################################

vflayoutmw_singleField2 <- function( vf, pwidth = 8.27,
                                    pheight = 11.69, margin = 0.25,
                                    filename = NULL,
                                    ffamily = "serif", sizetxt = 8,
                                    sizetxtSmall = 8, ffamilyvf = "serif", pointsize = 6,
                                    outerSymbol = "square", outerInch = 0.2,
                                    innerSymbol = "square", innerInch = outerInch / 2,
                                    lengthLines = 0, thicknessLines = 0) {
  
  # get normative values
  texteval <- "vfenv$nv"
  nv       <- eval( parse( text = texteval ) )
  
  if( nrow( vf ) > 1 ) {
    stop("Error! vf cannot have more than 1 row")
  }
  
  # open window with A4 page
  if( is.null( filename ) ) {
    device <- options( "device" )
    if( .Platform$OS.type == "unix" ) {
      if( Sys.info()["sysname"] == "Darwin" ) {
        options( device = "quartz" )
        dev.new( width = pwidth, height = pheight, dpi = 85 )
      } else {
        options( device = "x11" )
        dev.new( width = pwidth, height = pheight )
      }
    } else{
      options( device = "windows" )
      dev.new( width = pwidth, height = pheight, rescale = "fixed" )
    }
    options( device = device )
  } else {
    pdf( width = pwidth, height = pheight, file = filename )
  }
  
  # define the margins
  mwidth  <- pwidth  - 2 * margin
  mheight <- pheight - 2 * margin
  
  # create the layout of the printout
  printout <- createviewport( "printout", left = margin, top = margin, height = mheight, width = mwidth )

  ######################################################
  # first plot all graphs
  ######################################################
  texteval <- paste( vf$tperimetry, "locmap$", vf$tpattern, sep = "" )
  locmap   <- eval( parse( text = texteval ) )
  xrange <- max( locmap$xod ) - min( locmap$xod )
  yrange <- max( locmap$yod ) - min( locmap$yod )
  xmin   <- min( locmap$xod ) - 0.025 * xrange
  xmax   <- max( locmap$xod ) + 0.025 * xrange
  ymin   <- min( locmap$yod ) - 0.025 * yrange
  ymax   <- max( locmap$yod ) + 0.025 * yrange
  
  # total-deviation plot
  opar <- par( no.readonly = TRUE )
  par( fig = c( 0.35, 1, 0.666, 0.975 ) )
  vfplot( vf, plotType = "td", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, pointsize = pointsize)

    # raw sensitivities
  opar <- par( new = TRUE )
  par( fig = c( 0.35, 1, 0.3405, 0.6505 ) )
  vfplot( vf, plotType = "vf", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, pointsize = pointsize)

  # pattern-deviation plot
  opar <- par( new = TRUE )
  par( fig = c( 0.35, 1, 0.015, 0.325 ) )
  vfplot( vf, plotType = "pd", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, pointsize = pointsize)

  # color-code map
  par( new = TRUE )
  par( fig = c( 0.01, 0.21, 0.01, 0.16 ) )
  if( vfenv$nv$nvname == "nvsapmwcps" ) {
    colormapgraph( ncol = 5, notSeenAsBlack = FALSE, pointsize = 6, outerInch = 0.125 )
  } else {
    colormapgraph( ncol = 6, pointsize = pointsize, outerInch = 0.18 )
  }
  par( opar )

  ######################################################
  #
  #create the text elements in the printouts
  ######################################################
  # main info
  mainInfo         <- createviewport( "mainInfo",   left =  0.00, top =  0.00, width = 4.75, height = 0.40, pheight = mheight, pwidth = mwidth )
  # text for TD, PD, and sensitivity maps
  tdtext           <- createviewport( "tdtext",     left =  8.00, top =  0.00, width = 1.83, height = 0.22, pheight = mheight, pwidth = mwidth )
  senstext         <- createviewport( "senstext",   left =  8.00, top =  3.75, width = 1.83, height = 0.22, pheight = mheight, pwidth = mwidth )
  pdtext           <- createviewport( "pdtext",     left =  8.00, top =  7.50, width = 1.83, height = 0.22, pheight = mheight, pwidth = mwidth )
  # information about normative values and software version
  normvaltxt       <- createviewport( "normvaltxt", left =  6.37, top = 10.89, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  # test information
  typetxt      <- createviewport( "typetxt",     left =  0.05, top =  0.50, width = 1.40, height = 0.40, pheight = mheight, pwidth = mwidth )
  infotxt      <- createviewport( "infotxt",     left =  0.05, top =  0.75, width = 2.95, height = 0.40, pheight = mheight, pwidth = mwidth )
  labelreliab  <- createviewport( "labelreliab", left =  0.05, top =  1.10, width = 1.20, height = 0.65, pheight = mheight, pwidth = mwidth )
  
  reliability  <- createviewport( "reliability", left =  0.90, top =  1.10, width = 0.60, height = 0.65, pheight = mheight, pwidth = mwidth )
  label        <- createviewport( "label",       left =  0.05, top =  1.65, width = 0.50, height = 1.10, pheight = mheight, pwidth = mwidth )
  global       <- createviewport( "global",      left =  0.20, top =  1.65, width = 0.70, height = 1.10, pheight = mheight, pwidth = mwidth )
  pval         <- createviewport( "pval",        left =  1.00, top =  1.65, width = 1.00, height = 1.10, pheight = mheight, pwidth = mwidth )
  labeladd     <- createviewport( "labeladd",    left =  0.05, top =  2.50, width = 1.00, height = 1.10, pheight = mheight, pwidth = mwidth )
  addinfo      <- createviewport( "addinfo",     left =  1.50, top =  2.50, width = 1.00, height = 1.10, pheight = mheight, pwidth = mwidth )
  labelcomments <- createviewport( "labelcomments", left =  0.05, top =  3.40, width = 1.00, height = 1.10, pheight = mheight, pwidth = mwidth ) 
  comments     <- createviewport( "comments",    left =  0.05, top =  3.55, width = 1.00, height = 1.10, pheight = mheight, pwidth = mwidth )
  
  # create the list and then generate the tree and "push" it
  list <- vpList( mainInfo, tdtext, senstext, pdtext, normvaltxt,
                  typetxt, infotxt, labelreliab, reliability, label, global, pval, labeladd, addinfo, labelcomments, comments )
  tree <- vpTree( printout, list )
  
  pushViewport( tree )
  # map info
  seekViewport( "tdtext" )
  grid.text( "Total\nDeviation", x = 0.0, y = 1.0, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = 1.5 * sizetxt, fontface = "bold" ) )
  seekViewport( "senstext" )
  grid.text( "Sensitivity", x = 0.0, y = 1.0, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = 1.5 * sizetxt, fontface = "bold" ) )
  seekViewport( "pdtext" )
  grid.text( "Pattern\nDeviation", x = 0.0, y = 1.0, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = 1.5 * sizetxt, fontface = "bold" ) )
  
  ######################################################
  # perimetry information
  ######################################################
  seekViewport( "mainInfo" )
  text <- "Static Automated Perimetry. Single field analysis."
  # ID
  text <- paste( text, "Subject ID: ", sep = "\n" )
  text <- paste( text, vf$id, ",", sep = "" )
  # age
  text <- paste( text, " age: ", round( vf$sage ), ",", sep = "" )
  # eye
  texteye <- paste( "eye:", vf$seye, sep = " " )
  if( vf$seye == "OD" ) {
    texteye <- paste( texteye, "(right)", sep = " " )
  } else if ( vf$seye == "OS" ) {
    texteye <- paste( texteye, "(left)", sep = " " )
  } else {
    texteye <- paste( texteye, "(which?)", sep = " " )
  }
  text <- paste( text, texteye, sep = " " )
  grid.text( text, x = 0.0, y = 1.0, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = 1.5 * sizetxt, fontface = "bold" ) )
  
  ######################################################
  # Details about printouts
  ######################################################
  # algorithm
  textalgorithm <- "ZEST"
  
  seekViewport( "typetxt" )
  
  if( vf$tpattern == "p24d2" ) {
    textpattern <- "Central 24-2, Size III"
  } else if( vf$tpattern == "p30d2" ) {
    textpattern <- "Central 30-2, Size III"
  } else if( vf$tpattern == "p10d2" ) {
    textpattern <- "Central 10-2, Size III"
  } else if( vf$tpattern == "p24d2v" ) {
    textpattern <- "Central 24-2, Size V"
  } else if( vf$tpattern == "p30d2v" ) {
    textpattern <- "Central 30-2, Size V"
  } else if( vf$tpattern == "p10d2v" ) {
    textpattern <- "Central 10-2, Size V"
  } else if( vf$tpattern == "p30d1" ) {
    textpattern <- "Central 30-1, Size III"
  } else if( vf$tpattern == "p30d1v" ) {
    textpattern <- "Central 30-1, Size V"
  } else if (vf$tpattern == "peripheral") {
    textpattern <- "Peripheral, Size III"
  } else if (vf$tpattern == "peripheralv") {
    textpattern <- "Peripheral, Size V"
  } else if (vf$tpattern == "peripheralvi") {
    textpattern <- "Peripheral, Size VI"
  } else if( vf$tpattern == "pPC26" ) {
    textpattern <- "Polar Central 26, Size III"
  } else if( vf$tpattern == "pPC26v" ) {
    textpattern <- "Polar Central 26, Size V"
  } else if( vf$tpattern == "pPC26vi" ) {
    textpattern <- "Polar Central 26, Size VI"
  } else if( vf$tpattern == "pPeri_v2" ) {
    textpattern <- "Polar Peripheral, Size III"
  } else if(any(vf$tpattern == c("pPeriv","pPeriv_v2"))) {
    textpattern <- "Polar Peripheral, Size V"
  } else if( vf$tpattern == "pPerivi_v2" ) {
    textpattern <- "Polar Peripheral, Size VI"
  } else if( vf$tpattern == "pPT" ) {
    textpattern <- "Polar Total, Size III"
  } else if( vf$tpattern == "pPTv" ) {
    textpattern <- "Polar Total, Size V"
  } else if( vf$tpattern == "pPTvi" ) {
    textpattern <- "Polar Total, Size VI"
  } else {
    textpattern <- "Unknown"
  }
  
  text <- paste( textpattern, textalgorithm, sep = ", Algorithm ")
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ),
             gp = gpar( fontfamily = ffamily, fontsize = 1.2 * sizetxt, fontface = "bold") )
  
  text <- paste( textpattern, textalgorithm, sep = ", Algorithm ")
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ),
             gp = gpar( fontfamily = ffamily, fontsize = 1.2 * sizetxt , fontface = "bold") )
  
  ######################################################
  # Details about printouts
  ######################################################
  seekViewport( "normvaltxt" )
  
  text <- paste( "norm vals: ", visualFields::vfenv$nv$nvname, sep = "" )
  text <- paste( text, substr( packageDescription( "visualFields" )$Date, 1, 4 ), sep = "\n" )
  text <- paste( text, "visualFields", packageDescription( "visualFields" )$Version, sep = " " )
  grid.text( text, x = 1.00, y = 0.00, just = c( "right", "bottom" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxtSmall ) )
  
  ######################################################
  # visual-field test results
  ######################################################
  # subject and test information test
  seekViewport( "infotxt" )
  
  timetxt <- substr( vf$ttime, 1, 5 )
  if( substr( timetxt, 1, 1 ) == "0" ) substr( timetxt, 1, 1 ) <- ""
  text <- paste( "Date:", format( vf$tdate, "%m/%d/%Y" ), "at", timetxt, sep = " " )
  # duration and pause of test
  timetxt         <- substr( vf$sduration, 4, nchar( vf$sduration ) )
  if( timetxt != "59:59" ) {
    if( substr( timetxt, 1, 1 ) == "0" ) substr( timetxt, 1, 1 ) <- ""
    text <- paste( text, paste( "Duration: ", timetxt, sep = " " ), sep = "\n" )
  }
  timetxt         <- substr( vf$spause, 4, nchar( vf$sduration ) )
  if( timetxt != "59:59" ) {
    if( substr( timetxt, 1, 1 ) == "0" ) substr( timetxt, 1, 1 ) <- ""
    text <- paste( text, paste( ", pause: ", timetxt, sep = "" ), sep = "" )
  }
  text <- paste( text, "", sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  # false positives, negatives, and fixation losses
  seekViewport( "labelreliab" )
  
  text <- "FP"
  text <- paste( text, "FN", sep = "\n" )
  text <- paste( text, "FL", sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  seekViewport( "reliability" )
  
  sfp <- paste( sprintf( "%.1f", round( 1000 * vf$sfp ) / 10 ), "%", sep = " " )
  sfn <- paste( sprintf( "%.1f", round( 1000 * vf$sfn ) / 10 ), "%", sep = " " )
  sfl <- paste( sprintf( "%.1f", round( 1000 * vf$sfl ) / 10 ), "%", sep = " " )
  
  text <- sfp
  text <- paste( text, sfn, sep = "\n" )
  text <- paste( text, sfl, sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) ) 
  
  # global indices
  vfs  <- vfstats( vf )
  vfsp <- vfstatspmap( vfs )
  vfi  <- vfindex( vf, vfiset = vfidefault)
  vfip <- vfindexpmap( vfi )
  # general-height difference, if the used normative values have one.
  texteval <- paste( "vfenv$nv$", vf$tpattern, "_", vf$talgorithm, "$nvtdrank$mtdr", sep = "" )
  tdr <- NULL
  tdr <- eval( parse( text = texteval ) )
  if( !is.null( tdr ) ) {
    gh <- ghpostd( tdval( vf ) )
    gh <- paste( sprintf( "%.1f", round( 10 * gh ) / 10 ), "dB", sep = " " )
  }
  
  ms  <- paste( sprintf( "%.1f", round( 10 * vfs$msens ) / 10 ), "dB", sep = " " )
  md  <- paste( sprintf( "%.1f", round( 10 * vfs$mtdev ) / 10 ), "dB", sep = " " )
  psd <- paste( sprintf( "%.1f", round( 10 * vfs$spdev ) / 10 ), "dB", sep = " " )
  vfi <- paste( sprintf( "%.1f", round( 10 * vfi$mvfi  ) / 10 ), " %", sep = " " )
  
  seekViewport( "label" )
  
  text <- "MS"
  text <- paste( text, "MD", sep = "\n" )
  text <- paste( text, "PSD", sep = "\n" )
  text <- paste( text, "VFI", sep = "\n" )
  if( !is.null( tdr ) ) {
    text <- paste( text, "GH", sep = "\n" )
  }
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  seekViewport( "global" )
  
  text <- ms
  text <- paste( text, md, sep = "\n" )
  text <- paste( text, psd, sep = "\n" )
  text <- paste( text, vfi, sep = "\n" )
  if( !is.null( tdr ) ) {
    text <- paste( text, gh, sep = "\n" )
  }
  grid.text( text, x = 1.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  seekViewport( "pval" )
  
  text <- ""
  textp <- paste( "(p < ", vfsp$mtdev, " %)", sep = "" )
  text <- paste( text, textp, sep = "\n" )
  textp <- paste( "(p < ", vfsp$mtdev, " %)", sep = "" )
  text <- paste( text, textp, sep = "\n" )
  textp <- paste( "(p < ", vfip$mvfi, " %)", sep = "" )
  text <- paste( text, textp, sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  seekViewport( "labeladd" )
  text <- "Foveal Threshold:"
  text <- paste( text, "Total Presentations:", sep = "\n" )
  text <- paste( text, "Outlier Presentations:", sep = "\n" )
  text <- paste( text, "Mistakes:", sep = "\n" )
  text <- paste( text, "Deleted:", sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  seekViewport( "addinfo" )
  text <- vf$fovth
  text <- paste( text, vf$np, sep = "\n" )
  text <- paste( text, vf$outnp, sep = "\n" )
  text <- paste( text, vf$mistakes, sep = "\n" )
  text <- paste( text, vf$deleted, sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  seekViewport ("labelcomments")
  text <- "Comments:"
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  seekViewport( "comments" )
  text <- as.character( vf$comments )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
  # only if in save mode, then set device to off
  if( !is.null( filename ) ) {
    dev.off()
  }
  
}