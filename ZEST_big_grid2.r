# Perform a ZEST procedure at each location in a 24-2 pattern (right eye)
# using the HFA growth pattern to set the priors at each location.
#
# The prior pdf at each location is bimodal, with the guess from the HFA growth
# pattern as the mode of the normal part of the pdf, and the damaged part
# is fixed. The guess for the primary points (+-9,+-9) is passed in as a
# parameter to runBiModalZEST().
#
# The domain of the prior is -5:45 dB, which you might want to change.
#
# See example usage at end of file. 
#
# Author: Andrew Turpin & Luke Chong
# Date: Wed 12 Jun 2013 13:06:22 EST
#
# Modified 1 Feb 2015: made 18*30 matrix of locations (aturpin)
#                    : added 3 choices of grid
#         22 Feb 2015: added window for entering patient details and initiating test (lxchong)
#         23 Feb 2015: added FP & FN catch trials
#         24 Feb 2015: added test duration calculation (lxchong)
#                    : altered growth pattern for "Peripheral" grid
#                    : added trial-by-trial data (location, stim, response, response time)
#                      gets written to a .txt file
#         25 Feb 2015: added x and y to state for use by procedureWithGrowthPattern()
#          9 Mar 2015: reduced response window to 1000 ms
#         10 Mar 2015: added nextStimState parameter - tells projector where to position itself for next stimulus. 
#         11 Mar 2015: added pause and resume button.
#         16 Mar 2015: made adaptive response window (respWin)
#         18 Mar 2015: added new rule for stim choice so that a location can't be tested
#                      twice in a row except if it is the last location in a wave. 
#         08 Apr 2015: created new pauseAtStart button so test commences with a pause.
#                      operator needs to press button to initiate test.
#         10 Apr 2015: added drop down box in menu interface to select stim size.
#         12 Apr 2015: added  live test status window 
#         14 Apr 2015: added testing of outliers
#         27 May 2015: 30-1 grid added
#         12 Jun 2015: practice test added
#         29 Jul 2015: FP and FN catch trials removed from outlier test
#         18 Aug 2015: output files are saved in separate folders based on diagnosis
#         20 Aug 2015: Added new button when test paused, allowing termination of test without saving.
#         22 Aug 2015: Added new function where you can pause and delete the most recent presentations
#          1 Jan 2016: customized primary start value of BM priors based on MW normal data
#          2 Jan 2016: Front loaded most catch trials to the first minute of test.
#          6 Jan 2016: Added response time reliability measures to the live perimetrist display. 
#         13 Jan 2016: Updated input menu to include Rx, ORx and VA.
#         15 Jan 2016: Added ability to load details of previous patients at the input screen.
#         21 Jan 2016: Projector moves to next stimulus position whilst observer is responding (peripheral test only)
#         25 Jan 2016: Created adaptive interStimInterval based on observer's response times.
#         25 Apr 2016: Integrated procedure with the visualFields package
#         28 Apr 2016: Threshold Censoring and addition of new polar grids
#         22 July 2016: Added function to combine central and peripheral polar grids.
#         11 August 2016: Added Full From Prior option. Seeds unimodal prior based on previous test result.
#         13 August 2016: Changed outliers calculation rule
#         17 August 2016: Customized prior seeding for polar grids based on normative data collected.
#         25 August 2016: Added P-Central26 and P-Peripheral visualFields printouts
#           02 June 2018: Added new voronoi tessellation visualFields printout
#           02 June 2018: removed 16 locations from pPeriv test pattern
#           12 June 2018: added normative database for P-Total and added automatic generation of printout

rm(list=ls())
source("grids2.r")
source("growthPattern2.r")
source("pauseButton.r")
source("testStatusOutput.r")

# IMF: added parameters for initialization
bigWheel = FALSE;
resp_buzzer = 3;

#########################################################################
# INPUTS:
#   eye               - "right" or "left"
#   primaryStartValue - the mode of the normal part of the prior pdf for
#                       locations (+-9,+-9).
#   gridType          - "24-2", "30-2" or "Peripheral"
#   outlierFreq &     - If there are an outlierValFreq number of neighbouring locations where difference 
#   outlierValue        is >= outlierValue, that location is flagged as an outlier and is re-tested
#
#   If running in simulation mode you need to specify 
#     tt - an 109*181 matrix 24-2 field OD (blind spot on right) 
#     fpv - false positive rate in range [0,1]
#     fnv - false negative rate in range [0,1]
#   
# RETURNS: a list with
#    np = matrix of number of presenations for each location
#    ae = matrix of absolute error for each location
#    th = matrix of measured thresholds for each location
#    trues = matrix of input thresholds
#
# ALGORITHMS: Makes use of procedureWithGrowthPattern(...).
#
# 
#########################################################################
Zest242 <- function(eye="right", primaryStartValue=30, gridType="24-2",
                    minInterStimInterval=300,
                    tt=NA, fpv=0.00, fnv=0.00,outlierFreq=1,outlierValue=5,moveProjector = TRUE, retest=FALSE) {
    ####################################################################
    # Each location derives its start value from the average of all of the
    # immediate 9 neighbours that are lower than it.
    # Numbers should start at 1 and increase, not skipping any.
    ####################################################################
    
    if (gridType == "Peripheral") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.peripheral.coords)
            growthNext <- growth_lookup(grid.peripheral.coords)
        } else {
            growthPattern <- grid.flip(testPattern(grid.peripheral.coords))
            growthNext <- growth_lookup(grid.peripheral.coords, invert = TRUE)
        }
    } else if (gridType == "30-2") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.30.2.coords)
            growthNext <- growth_lookup(grid.30.2.coords)
        } else {
            growthPattern <- grid.flip(testPattern(grid.30.2.coords))
            growthNext <- growth_lookup(grid.30.2.coords, invert = TRUE)
        }
    } else if (gridType == "24-2") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.24.2.coords)
            growthNext <- growth_lookup(grid.24.2.coords)
        } else {
            growthPattern <- grid.flip(testPattern(grid.24.2.coords))
            growthNext <- growth_lookup(grid.24.2.coords, invert = TRUE)
        }
    } else if (gridType == "30-1") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.30.1.coords)
            growthNext <- growth_lookup(grid.30.1.coords)
        } else {
            growthPattern <- grid.flip(testPattern(grid.30.1.coords))
            growthNext <- growth_lookup(grid.30.1.coords, invert =  TRUE)
        }
    } else if (gridType == "G1") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.G1.coords)
            growthNext <- growth_lookup(grid.G1.coords)
        } else {
            growthPattern <- grid.flip(testPattern(grid.G1.coords))
            growthNext <- growth_lookup(grid.G1.coords, invert =  TRUE)
        }      
    } else if (gridType == "P-Total") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.PTotal.coords)
            growthNext <- growth_lookup(grid.PTotal.coords, searchExt = 12)
        } else {
            growthPattern <- grid.flip(testPattern(grid.PTotal.coords))
            growthNext <- growth_lookup(grid.PTotal.coords, searchExt = 12, invert =  TRUE)
        } 
    } else if (gridType == "P-Central10") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.PCentral10.coords)
            growthNext <- growth_lookup(grid.PCentral10.coords, searchExt = 5)
        } else {
            growthPattern <- grid.flip(testPattern(grid.PCentral10.coords))
            growthNext <- growth_lookup(grid.PCentral10.coords, searchExt = 5, invert =  TRUE)
        } 
    } else if (gridType == "P-Central26") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.PCentral26.coords)
            growthNext <- growth_lookup(grid.PCentral26.coords, searchExt = 12)
        } else {
            growthPattern <- grid.flip(testPattern(grid.PCentral26.coords))
            growthNext <- growth_lookup(grid.PCentral26.coords, searchExt = 12, invert =  TRUE)
        } 
    } else if (gridType == "P-Peripheral") {
        if (eye == "right") {
            growthPattern <- testPattern(grid.PPeri.coords)
            growthNext <- growth_lookup(grid.PPeri.coords, searchExt = 12)
        } else {
            growthPattern <- grid.flip(testPattern(grid.PPeri.coords))
            growthNext <- growth_lookup(grid.PPeri.coords, searchExt = 12, invert =  TRUE)
        } 
    } else if (gridType == "P-Edge") {
      if (eye == "right") {
        growthPattern <- testPattern(grid.PEdge.coords)
        growthNext <- growth_lookup(grid.PEdge.coords, searchExt = 16)
      } else {
        growthPattern <- grid.flip(testPattern(grid.PEdge.coords))
        growthNext <- growth_lookup(grid.PEdge.coords, searchExt = 16, invert =  TRUE)
      } 
    } else if (gridType == "practice") {
        growthNext <- NULL
        if (eye == "right") {
            growthPattern <- testPattern(grid.practice.coords)
        } else {
            growthPattern <- grid.flip(testPattern(grid.practice.coords))
        }
    } else {
        stop("Invalid gridType for Zest242")
    }

    onePriors <- apply(growthPattern, 1:2, function(x) ifelse(is.na(x), NA, ifelse(x == 1, primaryStartValue, NA)))

        # check tt is in the right format
    if (is.matrix(tt)) {
        if (ncol(tt) != 181 || nrow(tt) != 109) {
          stop("Zest242(): tt is not a 109*181 matrix")
        }
        z1 <- !is.na(growthPattern)
        z2 <- !is.na(tt)
        if (!all((z1 & z2) == z1))
            stop(paste("Zest242(): tt matrix is not the correct style for a",eye,"eye"))
    }


    ###############################
    # Set up fp and fn matrices
    ###############################
    z <- !is.na(growthPattern)
    fp <- z * matrix(fpv, nrow(growthPattern), ncol(growthPattern))
    fn <- z * matrix(fnv, nrow(growthPattern), ncol(growthPattern))
    
    
    ###############################################################
    # Seed prior from previous tests if retest == TRUE
    # prevTh - matrix of thresholds from previous test
    #
    ###############################################################
    
    if (retest == TRUE) {
      retestFile <- paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$dx,"_",details$gridType,
                         "_Grid_Size_",details$stimSizeRoman,"_vfPackage.csv")
    
      if (file.exists(retestFile)) {
        readRetest <- read.csv(retestFile)
      
        #Extract rows for previous test of patient for a given eye.
        #Check if previous test exists
        getPx <- intersect(grep(details$name,readRetest[,"id"]),grep(ifelse(details$eye == "left","OS","OD"),readRetest[,"seye"]))
        if (length(getPx) == 0) {
          stop("Unable to seed from prior. Previous test does not exist")
        }
      
        readRetest <- readRetest[getPx,]
        readRetest <- readRetest[nrow(readRetest),] #Extract most recent test
      
        prevTh <- matrix(NA,nrow = nrow(growthPattern),ncol = ncol(growthPattern),byrow=TRUE)
        retestIndex <- which(!is.na(growthPattern),arr.ind=TRUE) #find indices of test locations
        retestIndex <- retestIndex[order(retestIndex[,1]),] #order by row
        threshold <- readRetest[grep("L",names(readRetest))]
        if (details$gridType == "30-2") {threshold <- threshold[-which(threshold == 0)]} #remove blind spot locations for 30-2
      
        for (i in 1:length(threshold)) {
          prevTh[retestIndex[i,1],retestIndex[i,2]] <- as.numeric(threshold[i])
        }
      } else {
        stop("Unable to seed from prior. Previous test does not exist")
      }
    }
    
    #####################################################
    # Given guess, create a state for location (rw, cl)
    # where (rw,cl) is index into grid.* matrix
    # This version assumes
    #  1) bi-modal prior, constructed as in 
    #     Turpin et al  IOVS 44(11), November 2003. Pages 4787-4795.
    #  2) a domain of -5:45 dB
    #####################################################
    startF <- function(guess, rw, cl) {

        # censor thresholds
        if (any(details$gridType == c("P-Total","P-Central10","P-Central26","P-Peripheral","P-Edge"))) {
          minDomain <- 10
        } else {
          minDomain <- -5
        }
      
        domain <- minDomain:45 # note min(domain) should be <= 0 for bimodal pdf

        glaucomaPDF <- rep(0.001,length(domain))
        glaucomaPDF[1:10] <- c(rep(0.001,4),0.2, 0.3, 0.2, 0.15, 0.1, 0.02) #BM PRIOR
        #glaucomaPDF[1] <- 0.2 #custom
        #glaucomaPDF[1:6] <- c(0.2, 0.3, 0.2, 0.15, 0.1, 0.02) # BM2 PRIOR
        #uniformPDF <- rep(0.001,length(domain))
        #gaussPDF <- dnorm(domain,guess,6)
        healthyPDF  <- function (normalModePoint) {
            temp <- c(rep(0.001,50), 0.001, 0.03, 0.05, 0.1, 0.2, 0.3, 0.2, 0.05, 0.025, 0.01,rep(0.001,50))
            mode <- which.max(temp)
            return(temp[(mode-normalModePoint+domain[1]):(mode-normalModePoint+domain[length(domain)])])
        }
        makeBiModalPDF <- function(normalModePoint, weight, pdf.floor) {
          npdf <- healthyPDF(normalModePoint)
          cpdf <- npdf * weight + glaucomaPDF
          cpdf[which(cpdf < pdf.floor)] = pdf.floor 
          return (cpdf)
        }	
        
        if (retest == TRUE && details$fovea == FALSE) {
          glaucomaPDF <- rep(0.001,length(domain)) 
          prior <- makeBiModalPDF(max(prevTh[rw,cl], minDomain + 5),4,0.001)  #if prevTh is below censored th, set peak at censored th. 
        } else {
          prior <- makeBiModalPDF(round(guess),4,0.001)
        }
        
        #prior <- uniformPDF
        #prior <- gaussPDF
        #prior[which(prior < 0.001)] <- 0.001
        prior <- prior/sum(prior)
      
        ###### helps create the makeStim function needed for the OPI
        makeStimHelper <- function(x, y) {  # returns a function of (db,n)
            ff <- function(db, n) db+n
            body(ff) <- substitute(
                {s <- list(x=x, y=y, level=dbTocd(db,4000/pi),
                        size=as.numeric(details$stimSize),
                        color=.OpiEnv$O900$STIM_WHITE,
                        duration=200, responseWindow=1000)
                 
                 class(s) <- "opiStaticStimulus"
                return(s)
                }
                , list(x=x,y=y))
            return(ff)
        }
          
        ms <- makeStimHelper(-6-14*6 + (cl-1), 6+8*6 - (rw-1))
        state <- ZEST.start(domain=domain, prior=prior, makeStim=ms, 
                              minStimulus=max(0,minDomain),
                              maxStimulus=40,
                              tt=ifelse(!is.null(dim(tt)), tt[rw,cl], NA),
                              fpr=ifelse(!is.null(dim(tt)), fpv, NA),
                              fnr=ifelse(!is.null(dim(tt)), fnv, NA)
                              )
        state$x <- -6-14*6 + (cl-1)
        state$y <- 6+8*6 - (rw-1)
        
        return(state)
    }

    #####################################################
    # Given a state, step the procedure and return new state
    #####################################################
    stepF <- function(state, nextStimState=NULL) {
            #print(state$pdf)
            if (is.null(nextStimState))
                return(ZEST.step(state)$state)
            else
                return(ZEST.step(state, nextStim=nextStimState$makeStim(0,0))$state)
    }

    #####################################################
    # Given a state, return TRUE for finished, FALSE otherwise
    #####################################################
    stopF <- ZEST.stop

    #####################################################
    # Given a state, return c(threshold, num presentations)
    # Round -ve thresholds to 0
    #####################################################
    finalF <- function(state) {
        t <- ZEST.final(state)
        #print(state$stimuli)
        return(c(t, state$numPresentations))
    }
    
    ####################################################
    # function for foveal threshold testing
    ###################################################
    while (details$fovea == TRUE) {
      resFovea <- function () {
        state.fovea <- startF(primaryStartValue,55,91)
        state.fovea$opiParams$tt <- 30
        while (stopF(state.fovea) == FALSE) {
          Sys.sleep(runif(1, min=400, max=800)/1000)
          state.fovea <- stepF(state.fovea)         
        }
        return(finalF(state.fovea)[1])   
      }
      opiSetBackground(fixation=.OpiEnv$O900$FIX_CROSS, fixIntensity=30, lum=.OpiEnv$O900$BG_10)
      pauseAtStartFovea() 
      fovealTH <- resFovea()
      fovealTestComplete(fovealTH)
    }
    
    windows(700,250)
    opiSetBackground(fixation=.OpiEnv$O900$FIX_CENTRE, fixIntensity=30, lum=.OpiEnv$O900$BG_10)
    pauseAtStart()
    commence <<- Sys.time()
    res1 <- procedureWithGrowthPattern(details$startTime,growthPattern, growthNext, onePriors, startF, stepF, stopF, finalF,
            gridPat=growthPattern,
            respWinBuffer=250,
            FPLevel=dbTocd(55, 4000/pi), 
            FNDelta=10,
            FNPause=300,
            FNLocationThreshold=20,
            FPSize=as.numeric(details$stimSize),
            FNSize=as.numeric(details$stimSize),
            moveProj = moveProjector,
            minInterStimInt = minInterStimInterval)

    #z <- res1$t < 0
    tz <- res1$t
    #tz[z] <- 0
    
    #########################################################################
    # Calculate Difference between neighbours within each of the 4 quadrants
    # and identify outliers
    #########################################################################
    NeighbourDifference <- function(vf,outlierValue,search_extent = 10) {
      result <- matrix(NA,nrow(vf),ncol(vf))
      for (i in list(1:(0.5*nrow(vf)),(0.5*nrow(vf)+1):nrow(vf))) {
        for (j in list(((0.5*ncol(vf)+1):ncol(vf)),(1:(0.5*ncol(vf))))) {
          quadrant <- vf[i,j]
          quadResult <- matrix(NA,nrow(quadrant),ncol(quadrant))
          locIndex <- which(!is.na(quadrant),arr.ind=TRUE)
          for (rr in 1:nrow(locIndex)) {
            deltas <- NULL
            distances <- apply(locIndex, 1, function (t) {
              round(sqrt((t[2] - locIndex[rr,2])^2 + (t[1] - locIndex[rr,1])^2))
            })
            neighbours <- locIndex[(distances > 0 & distances <= search_extent),]
            neighbours <- matrix(neighbours,length(unlist(neighbours)) / 2,2)
            
            if (nrow(neighbours) == 0) {
              deltas <- NULL
            } else {
              if (nrow(neighbours) == 1) {
                outlierFreq <- 1
              } else if (any(nrow(neighbours) == c(2:4))) {
                outlierFreq <- 2
              } else if (any(nrow(neighbours) == c(5:6))) {
                outlierFreq <- 3
              } else if (nrow(neighbours) == 7) {
                outlierFreq <- 4
              } else if (nrow(neighbours) >= 8) {
                outlierFreq <- 5
              }
              deltas <- round(apply(neighbours, 1, function (t) {
                abs(quadrant[locIndex[rr,1],locIndex[rr,2]] - quadrant[t[1],t[2]])
              }))
              
              if (sum(deltas >= outlierValue, na.rm=TRUE) >= outlierFreq) {quadResult[locIndex[rr,1],locIndex[rr,2]] <- 1}  
            } 
          }
          result[i,j] <- quadResult
        }
      }
      return(result)
    }
    
    ###########################################################
    # Re-test outliers
    ###########################################################
    if (any(gridType == c("P-Total","P-Peripheral"))) {
      outliers <- NeighbourDifference(tz,outlierValue,search_extent = 14)
    } else if (gridType == "P-Central10") {
      outliers <- NeighbourDifference(tz,outlierValue,search_extent = 5)
    } else if (gridType == "P-Edge") {
      outliers <- NeighbourDifference(tz,outlierValue,search_extent = 16)
    } else {
      outliers <- NeighbourDifference(tz,outlierValue,search_extent = 10)
    }
    
    if (sum(outliers,na.rm=TRUE) > 0) {
      outlierPriors <- apply(outliers, 1:2, function(x) ifelse(is.na(x), NA, ifelse(x == 1, primaryStartValue, NA)))
      windows(width=700,height=250,ypos=300)
    
      res2 <- procedureWithGrowthPattern(details$startTime,outliers, gn = NULL, outlierPriors, startF, stepF, stopF, finalF, 
            gridPat=growthPattern,
            respWinBuffer=250,
            catchTrialLoadFreq = 5000,
            catchTrialFreq=5000,
            FPLevel=dbTocd(55, 4000/pi), 
            FNDelta=10,
            FNPause=300,
            FNLocationThreshold=20,
            FPSize=as.numeric(details$stimSize),
            FNSize=as.numeric(details$stimSize),
            moveProj = moveProjector,
            minInterStimInt = minInterStimInterval)
      
      #q <- res2$t < 0
      tq <- res2$t
      #tq[q] <- 0
    } else {res2 <- list(thOutliers=matrix(NA,nrow = nrow(growthPattern),ncol = ncol(growthPattern),byrow=TRUE),
                         npOutliers=matrix(NA,nrow = nrow(growthPattern),ncol = ncol(growthPattern),byrow=TRUE))}
    
    if (exists("fovealTH")) {
      res <- list(np=res1$n, ae=abs(tz-tt), th=res1$t,thOutliers=res2$t,npOutliers=res2$n,thFovea=fovealTH,trues=tt,fp_counter=res1$fpc,fn_counter=res1$fnc,fp_outliers=res2$fpc,fn_outliers=res2$fnc,rt=res1$rt,rtOutliers=res2$rt)
    } else {
      res <- list(np=res1$n, ae=abs(tz-tt), th=res1$t,thOutliers=res2$t,npOutliers=res2$n,trues=tt,fp_counter=res1$fpc,fn_counter=res1$fnc,fp_outliers=res2$fpc,fn_outliers=res2$fnc,rt=res1$rt,rtOutliers=res2$rt)
    }
    return(res)
}

###########################################################################################
# Function which combines two separate test results into one
###########################################################################################
combine <- function(test1,test2) {
  index <- which(!is.na(test1[["np"]])) #find the index of central grid test locations
  z <- vector("list",length(test1)) #create empty list
  names(z) <- names(test1) # assign list names

  for (om in c("np","ae","th","thOutliers","npOutliers","trues")) {
    test2[[om]][index] <- test1[[om]][index]  #add central grid to peripheral
    z[[om]] <- test2[[om]]
  }
  
  for (om2 in c("fp_counter","fn_counter","fp_outliers","fn_outliers","rt","rtOutliers")) {
    z[[om2]] <- c(test1[[om2]],test2[[om2]])
  }
  
  if (!is.null(z1$thFovea) || !is.null(z2$thFovea)) {
    z[["thFovea"]] <- c(z1$thFovea,z2$thFovea)
  }
  
  return(z)
}

###########################################################################################
# Function which sets the primary start value for BM priors
# depending on test pattern/target size combination.
# NOTE: Procedure breaks down if PSV > 30. PSV must be 30 dB or lower.
###########################################################################################
setPSV <- function (grid,size) {
  if (size == "III") {
      if (grid == "30-2") {PSV <- 28}
      else if (grid == "P-Central10") {PSV <- 30}
      else {PSV <- 30} 
  } else if (size == "V") {
      if (any(grid == c("30-2","30-1","P-Central26"))) {PSV <- 33}
      else if (grid == "Peripheral") {PSV <- 30}
      else if (grid == "P-Peripheral") {PSV <- 31}
      else {PSV <- 32}
  } else if (size == "VI") {
      if (any(grid == c("Peripheral","P-Edge"))) {PSV <- 32}
      else {PSV <- 32}
  }
  return(PSV)
}

######################################################################
#
# create data frame containing px details and write to .csv file
# INPUT: 
#
# filename        - name of the output file
#
######################################################################
writeFile <- function (filename = paste(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$name,"_",details$dx,"_",details$gridType,"_",details$stimSizeRoman,"_",details$eye,"Eye_",details$date,"_",details$startTime,".csv",sep="")) {
  px_info <- data.frame(stringsAsFactors=FALSE)
  labels <- c("#ID:","#Age:","#Diagnosis:","#Manifest Refraction:","#Over-refraction:","#VA:","#Comments:","#Grid Type:","#Stimulus Size:","#Eye:")
  px_info <- data.frame(cbind(c(px_info,labels)))
  for (i in 1:10) {
    px_info[i,2] <- paste(details[i],sep="")
  }
  cat(paste("#",date(),sep=""), "\n",  file=filename, append=TRUE)
  px_info <- data.frame(lapply(px_info, as.character), stringsAsFactors=FALSE)
  write.table(px_info,filename,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)

   ######################################################################
   #
   # append threshold and presentation number data to px details file
   #
   #########################################################################

  cat(paste("#Test Duration: ",format(.POSIXct(difftime(terminate,commence,units="secs") - timePaused),"%M:%S"),sep=""),"\n",file=filename,append=TRUE)
  cat(paste("#Total Presentations: ", sum(c(unlist(z$np),unlist(z$npOutliers)),na.rm=TRUE),sep=""),"\n",file=filename,append=TRUE)
  cat(paste("#FP errors: ",sum(z$fp_counter,na.rm=TRUE),"/",length(z$fp_counter)," (",signif(sum(z$fp_counter,na.rm=TRUE)/length(z$fp_counter)*100,digits=3),"%)",sep=""), "\n",file=filename,append=TRUE)
  cat(paste("#FN errors: ",sum(z$fn_counter,na.rm=TRUE),"/",length(z$fn_counter)," (",signif(sum(z$fn_counter,na.rm=TRUE)/length(z$fn_counter)*100,digits=3),"%)",sep=""), "\n",file=filename,append=TRUE)
  cat(paste("#Full From Prior: ",details$retest),"\n",file=filename,append=TRUE)
  
   if (!is.null(z$thFovea)) {
     cat(paste("#Foveal Threshold: ",round(z$thFovea)," dB", sep=""),"\n",file=filename,append=TRUE)
   } 

  cat("\n#Thresholds\n", file=filename, append=TRUE)
  write.table(round(z$th[apply(z$th,1,function (x) !all(is.na(x))),apply(z$th,2,function (x) !all(is.na(x)))]), file=filename, append=TRUE,row.names=FALSE,col.names=FALSE,sep=",")
  if (any(!is.na(z$thOutliers))) {
    cat("\n#Outlier Thresholds\n", file=filename, append=TRUE)
    
    write.table(round(z$thOutliers[apply(z$th,1,function (x) !all(is.na(x))),apply(z$th,2,function (x) !all(is.na(x)))]), file=filename, append=TRUE,row.names=FALSE,col.names=FALSE,sep=",")
   }
  cat("\n#Num Presentations\n", file=filename, append=TRUE)
  write.table(z$np[apply(z$np,1,function (x) !all(is.na(x))),apply(z$np,2,function (x) !all(is.na(x)))], file=filename, append=TRUE,row.names=FALSE,col.names=FALSE,sep=",")
  cat("\n", file=filename, append=TRUE)
}

##########################################################################################################################
#
# a second writeFile function which writes a .csv file with information displayed in columns
# 
#
##########################################################################################################################
writeFile2 <- function (details,filename = paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$dx,"_",details$gridType,"_Grid_Size_",details$stimSizeRoman,".csv")) {
  
  #################################################################################
  # Arrange thresholds into an appropriate format
  #
  ################################################################################
  thTable <- function (details) {
    t <- NULL
    if (details$eye == "left") {
      th <- grid.flip(z$th)
    } else {th <- z$th}
    for (rw in 1:nrow(th)) {
      for (cl in 1:ncol(th)) {
        if (!is.na(th[rw,cl])) {
          t <- cbind(t,round(th[rw,cl]))
          colnames(t)[ncol(t)] <- paste(-6-14*6 + (cl-1),6+8*6 - (rw-1),sep=",")
        }
      }
    }
    return(t)
  }
  ###############################################################################
  #
  # Compile the output spreadsheet
  #
  ###############################################################################
  output <- data.frame(ID=details$name,
                       Age = details$age,
                       Diagnosis = details$dx,
                       Manifest_Rx = details$MRx,
                       ORx = details$OR,
                       VA = details$VA,
                       Comments = details$comments,
                       Grid = paste0("'",details$gridType,"'"),
                       Stimulus_Size = details$stimSizeRoman, 
                       Eye = ifelse(details$eye == "right","RE","LE"),
                       Date = details$date,
                       Time = details$startTime,
                       Test_Duration = format(.POSIXct(difftime(terminate,commence,units="secs") - timePaused),"%M:%S"),
                       Presentations = sum(c(unlist(z$np),unlist(z$npOutliers)),na.rm=TRUE),
                       FP_seen = sum(z$fp_counter,na.rm=TRUE),
                       FP_total = length(z$fp_counter),
                       FN_missed = sum(z$fn_counter,na.rm=TRUE),
                       FN_total = length(z$fn_counter),
                       Foveal_Th = ifelse(!is.null(z$thFovea),round(z$thFovea),NA))
  
  output <- cbind(output,thTable(details))
  
  if (any(!is.na(z$thOutliers))) {
    if (details$eye == "left"){
      outs <- grid.flip(z$thOutliers)
    } else {outs <- z$thOutliers}
  
    outlierIndex <- which(!is.na(outs),arr.ind=TRUE)
    for (i in 1:nrow(outlierIndex)) {
      round(outs[outlierIndex[i,][1],outlierIndex[i,][2]])
      output[,paste(-6-14*6 + (outlierIndex[i,2]-1),6+8*6 - (outlierIndex[i,1]-1),sep=",")] <- paste(output[,paste(-6-14*6 + (outlierIndex[i,2]-1),6+8*6 - (outlierIndex[i,1]-1),sep=",")], "(",round(outs[outlierIndex[i,][1],outlierIndex[i,][2]]),")",sep="")
    }
  }
  
  if (file.exists(filename)) {
    write.table(output,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
  } else {
    write.table(output,file=filename,append=TRUE,row.names=FALSE,sep=",")
  }
}

###########################################################################################
# Function which writes a .csv file in a format that can be used in the "visualFields"
# package and creates printout using visualFields package
#
###########################################################################################
writeFile3 <- function (details,filename = paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$dx,"_",details$gridType,"_Grid_Size_",details$stimSizeRoman,"_vfPackage.csv")) {
  
  if (details$gridType == "24-2") {
    if (details$stimSizeRoman == "III") {pattern <- "p24d2"}
    if (details$stimSizeRoman == "V") {pattern <- "p24d2v"}
    if (details$stimSizeRoman == "VI") {pattern <- "p24d2vi"}
  }
  if (details$gridType == "30-2") {
    if (details$stimSizeRoman == "III") {pattern <- "p30d2"}
    if (details$stimSizeRoman == "V") {pattern <- "p30d2v"}
    if (details$stimSizeRoman == "VI") {pattern <- "p30d2vi"}
  }
  if (details$gridType == "30-1") {
    if (details$stimSizeRoman == "III") {pattern <- "p30d1"}
    if (details$stimSizeRoman == "V") {pattern <- "p30d1v"}
    if (details$stimSizeRoman == "VI") {pattern <- "p30d1vi"}
    }
  if (details$gridType == "Peripheral") {
    if (details$stimSizeRoman == "III") {pattern <- "peripheral"}
    if (details$stimSizeRoman == "V") {pattern <- "peripheralv"}
    if (details$stimSizeRoman == "VI") {pattern <- "peripheralvi"}
  }
  if (details$gridType == "G1") {
    if (details$stimSizeRoman == "III") {pattern <- "pG1"}
    if (details$stimSizeRoman == "V") {pattern <- "pG1v"}
    if (details$stimSizeRoman == "VI") {pattern <- "pG1vi"}
  }
  if (details$gridType == "P-Total") {
    if (details$stimSizeRoman == "III") {pattern <- "pPT"}
    if (details$stimSizeRoman == "V") {pattern <- "pPTv"}
    if (details$stimSizeRoman == "VI") {pattern <- "pPTvi"}
  } 
  if (details$gridType == "P-Central10") {
    if (details$stimSizeRoman == "III") {pattern <- "pPC10"}
    if (details$stimSizeRoman == "V") {pattern <- "pPC10v"}
    if (details$stimSizeRoman == "VI") {pattern <- "pPC10vi"}
  }
  if (details$gridType == "P-Central26") {
    if (details$stimSizeRoman == "III") {pattern <- "pPC26"}
    if (details$stimSizeRoman == "V") {pattern <- "pPC26v"}
    if (details$stimSizeRoman == "VI") {pattern <- "pPC26vi"}
  } 
  if (details$gridType == "P-Peripheral") {
    if (details$stimSizeRoman == "III") {pattern <- "pPeri_v2"}
    if (details$stimSizeRoman == "V") {pattern <- "pPeriv_v2"}
    if (details$stimSizeRoman == "VI") {pattern <- "pPerivi_v2"}
  }
  if (details$gridType == "P-Edge") {
    if (details$stimSizeRoman == "III") {pattern <- "pEdge"}
    if (details$stimSizeRoman == "V") {pattern <- "pEdgev"}
    if (details$stimSizeRoman == "VI") {pattern <- "pEdgevi"}
  }
  
    
  output <- data.frame(id=details$name,
                     tperimetry = "sap",
                     talgorithm = "zest",
                     tpattern = pattern,
                     tdate = format(Sys.Date(),"%m/%d/%Y"),
                     ttime = paste0(unlist(strsplit(details$startTime,".",fixed=TRUE))[1],":",unlist(strsplit(details$startTime,".",fixed=TRUE))[2],":",unlist(strsplit(details$startTime,".",fixed=TRUE))[3]),
                     stype = details$dx,
                     sage = details$age,
                     seye = ifelse(details$eye == "right","OD","OS"),
                     sbsx = 15,
                     sbsy = -1,
                     sfp = round(sum(z$fp_counter,na.rm=TRUE)/length(z$fp_counter),digits=2), 
                     sfn = round(sum(z$fn_counter,na.rm=TRUE)/length(z$fn_counter),digits=2),
                     sfl = 0,
                     sduration = format(.POSIXct(difftime(terminate,commence,units="secs") - timePaused,tz="GMT"),"%H:%M:%S"),
                     spause = format(.POSIXct(timePaused,tz="GMT"),"%H:%M:%S"),
                     fovth = ifelse(!is.null(z$thFovea),round(z$thFovea),NA),
                     np = sum(c(unlist(z$np),unlist(z$npOutliers)),na.rm=TRUE),
                     outnp = sum(unlist(z$npOutliers),na.rm=TRUE),
                     mistakes = mistakes,
                     deleted = deletes, 
                     comments = details$comments)

  if (details$eye == "left") {
      th <- grid.flip(z$th)
  } else {th <- z$th}
  thresholds <- matrix(round(t(th)[!is.na(t(th))]),1,sum(!is.na(th)))
  
  if (any(!is.na(z$thOutliers))) {
    if (details$eye == "left"){
      outs <- grid.flip(z$thOutliers)
    } else {outs <- z$thOutliers}
   
    outliers <- t(outs)[which(!is.na(t(th)))]
    for (i in 1:length(outliers)) {
      if (!is.na(outliers[i])) {
        thresholds[i] <- round(outliers[i])
      }
    }
  }
  
  if (details$gridType == "30-2") { # add blind spot locations to 30-2 grid
  thresholds <- c(thresholds[1:35],0,thresholds[36:44],0,thresholds[45:length(thresholds)])
  thresholds <- matrix(thresholds,1,length(thresholds))
  }
  colnames(thresholds) <- sapply(1:ncol(thresholds), function (x) {paste0("L",x)})
  finalOutput <- cbind(output,thresholds)
    
  if (file.exists(filename)) {
    write.table(finalOutput,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE,sep=",",na="")
  } else {
    write.table(finalOutput,file=filename,append=TRUE,row.names=FALSE,sep=",",na="")
  }

}

#####################################################################################################
# Compile a patient database in order to recall existing patient details for the input menu.
#
#
#####################################################################################################
px_database <- function (details) {
  output <- data.frame(ID=details$name,
                       Age = details$age,
                       Diagnosis = details$dx,
                       Manifest_Rx = details$MRx,
                       ORx = details$OR,
                       VA = details$VA,
                       Comments = details$comments,
                       Grid = details$gridType,
                       Stimulus_Size = details$stimSize, 
                       Eye = details$eye)
  
  if (file.exists("patient_list.txt")) {
    PL <- read.csv("patient_list.txt")
    if (any(grepl(paste0("^",details$name,"$"),PL$ID))) {
      if (!any(PL$Eye[which(PL$ID == details$name)] == details$eye)) {
        write.table(output,file="patient_list.txt",append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
      }
    } else {
        write.table(output,file="patient_list.txt",append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
    }
  } else {
    write.table(output,file="patient_list.txt",append=TRUE,row.names=FALSE,sep=",")
  }
}

####################################
# Example Usage
####################################
require(OPI)
chooseOpi("Octopus900")
source("query_patient_details.r")

    # extra opiInitialize to light up bowl before procedure starts
opiInitialize(eyeSuiteSettingsLocation="C:/ProgramData/Haag-Streit/EyeSuite/",eye="right",bigWheel=bigWheel,resp_buzzer = resp_buzzer, zero_dB_is_10000_asb = FALSE)
opiClose()

gRunning <- TRUE

details <- practiceQuery()

while (details$practice == TRUE) {
  gRunning <- TRUE
  opiInitialize(eyeSuiteSettingsLocation="C:/ProgramData/Haag-Streit/EyeSuite/",eye=details$eye,bigWheel=bigWheel,resp_buzzer = resp_buzzer, zero_dB_is_10000_asb = FALSE)
  Zest242(eye=details$eye, primaryStartValue=30, gridType="practice",outlierValue=15,outlierFreq=1)
  tkdestroy(tt)
  pracTestComplete()
  dev.off()
  details <- practiceQuery()
  opiClose()
}

gRunning <- TRUE # reset gRunning in case practice test was terminated early
details <- inputs()
if (dir.exists(paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman)) == FALSE) {dir.create(paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman),recursive = TRUE)}
opiInitialize(eyeSuiteSettingsLocation="C:/ProgramData/Haag-Streit/EyeSuite/",eye=details$eye,bigWheel=bigWheel, resp_buzzer = resp_buzzer, zero_dB_is_10000_asb = FALSE)
PSV <- setPSV(details$gridType,details$stimSizeRoman)

if (details$gridType == "P-Total") {
  z1 <- Zest242(eye=details$eye, primaryStartValue=PSV, gridType="P-Central26",outlierValue=7,minInterStimInterval=0,moveProjector = TRUE,retest=details$retest)
  tkdestroy(tt)
  graphics.off()
  details$fovea <- FALSE
  z2 <- Zest242(eye=details$eye, primaryStartValue=PSV, gridType="P-Peripheral",outlierValue=7,minInterStimInterval=0,moveProjector = TRUE,retest=details$retest)
  z <- combine(z1,z2)
} else {
  z <- Zest242(eye=details$eye, primaryStartValue=PSV, gridType=details$gridType,outlierValue=7,minInterStimInterval=0,moveProjector = TRUE,retest=details$retest)
}

terminate <- Sys.time()
tkdestroy(tt)  # closes the pause button upon completion of the test
opiClose()
graphics.off()

if (gRunning) {
  windows(900,350)
  testStatusFinal(z)
  testComplete()
}

if (gRunning) {
  pdf(file = paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$name,"_",details$dx,"_",details$gridType,"_",details$stimSizeRoman,"_",details$eye,"Eye_",details$date,"_",details$startTime,".pdf"),width=14,height=6)
  testStatusFinal(z)
  dev.off()
  comments <- finalComments()
  details$comments <- paste(details$comments,comments,sep=" ")
  writeFile()
  writeFile2(details)
  writeFile3(details)
  px_database(details)

    if (any(details$gridType == c("30-2","30-1","24-2","Peripheral","P-Peripheral","P-Central26","P-Total"))) {  
    ###################################################################################################
    # Create printout of data using visualFields package
    ###################################################################################################
    library(visualFields)
    ################################################################################
    # load patches to visualFields, as new normative values, locations map, etc here
    ################################################################################
    load("nvsapmwcps.rda")
    load( "vfsettingsmw.rda" )
    source( "vflayoutmw_singleField.r" )
    source( "vflayoutmw_singleField2.r" )

    #CARE!!! set appropriate normative values
    setnv( "nvsapmwcps" )

    #load data
    filename <- paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$dx,"_",details$gridType,"_Grid_Size_",details$stimSizeRoman,"_vfPackage.csv")
    loadfile <- read.csv(filename)
    vf <- loadvfcsv( filename = filename, patternMap = eval(parse(text = paste0("saplocmap$",as.character((tail(loadfile$tpattern,1)))))))

    #generate unique file name for printout
    fname <- paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$name,"_",details$dx,"_",details$gridType,"_",details$stimSizeRoman,"_",details$eye,"Eye_",details$date,"_",details$startTime,"_visualFields.pdf")
    fname2 <- paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$name,"_",details$dx,"_",details$gridType,"_",details$stimSizeRoman,"_",details$eye,"Eye_",details$date,"_",details$startTime,"_visualFields2.pdf")

    #save printout
    vflayoutmw_singleField(vf[nrow(vf),], filename = fname)
    vflayoutmw_singleField2(vf[nrow(vf),], filename = fname2)
   }
} else {
  file.remove(paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman,"/",details$name,"_",details$dx,"_",details$grid,"_",details$stimSizeRoman,"_",details$eye,"Eye_",details$date,"_",details$startTime,"_stimResponses.txt"))
}

opiClose()
#source('ZEST_big_grid_simulation.r')  # simulation test
#source("censored_PDF_simulation.r")