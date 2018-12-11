####################################
# for simulation
####################################
require(OPI)
chooseOpi("SimHenson")
#chooseOpi("SimGaussian")
#chooseOpi("SimYes")
#chooseOpi("SimNo")
#opiInitialize(sd=1)
opiInitialize(type="C", cap=6, display=NULL,maxStim = 4000/pi)
#opiInitialize(sd=1)
source("query_patient_details.r")
gRunning <- TRUE

Peripheral_TT <- testPattern(grid.peripheral.coords)
Peripheral_TT[14815] <- 30

TT.30.1 <- testPattern(grid.30.1.coords)
TT.30.1[which(TT.30.1 == 2)] <- 40

TT.30.2 <- testPattern(grid.30.2.coords)
TT.30.2[!is.na(TT.30.2)] <- 30
#TT.30.2[which(TT.30.2 == 1)] <- 35
#TT.30.2[which(TT.30.2 == 2)] <- 33
#TT.30.2[which(TT.30.2 == 3)] <- 31
#TT.30.2[which(TT.30.2 == 4)] <- 29

TT.24.2 <- testPattern(grid.24.2.coords)
TT.24.2[which(TT.24.2 == 1)] <- 30

TT.G1 <- testPattern(grid.G1.coords)
TT.G1[which(TT.G1 == 1)] <- 30

TT.PTotal <- testPattern(grid.PTotal.coords)
#TT.PTotal[which(TT.PTotal == 1)] <- 30
TT.PTotal[which(is.na(TT.PTotal) == FALSE)] <- 30

TT.PCentral10 <- testPattern(grid.PCentral10.coords)
#TT.PCentral10[which(TT.PCentral10 == 2)] <- 30
TT.PCentral10[which(is.na(TT.PCentral10) == FALSE)] <- 30

TT.PCentral26 <- testPattern(grid.PCentral26.coords)
#TT.PCentral26[which(TT.PCentral26 == 2)] <- 30
TT.PCentral26[which(is.na(TT.PCentral26) == FALSE)] <- 30
#TT.PCentral26[7027] <- 1

TT.PEdge <- testPattern(grid.PEdge.coords)
TT.PEdge[which(TT.PEdge == 2)] <- 30

TT.PPeri <- testPattern(grid.PPeri.coords)
#TT.PPeri[which(TT.PPeri == 2)] <- 30
TT.PPeri[which(is.na(TT.PPeri) == FALSE)] <- 30
TT.PPeri[5493] <- 1


tt.practice <- testPattern(grid.practice.coords)

###############################################################################################
#
# Begin test!!
#
##############################################################################################
details <- practiceQuery()

while (details$practice == TRUE) {
  Zest242(eye=details$eye, primaryStartValue=30, gridType="practice",outlierValue=8,outlierFreq=2,tt=tt.practice)
  tkdestroy(tt)
  pracTestComplete()
  dev.off()
  details <- practiceQuery()
}

details <- inputs()
if (dir.exists(paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman)) == FALSE) {dir.create(paste0(details$dx,"/",details$gridType," ",details$stimSizeRoman),recursive = TRUE)}
np <- NULL
PSV <- setPSV(details$gridType,details$stimSizeRoman)
#trials <- 500
#for ( i in 1:trials) {
#print(i)

if (details$gridType == "Peripheral") {
  TT <- Peripheral_TT
} else if (details$gridType == "30-1") {
  TT <- TT.30.1
} else if (details$gridType == "30-2") {
  TT <- TT.30.2
} else if (details$gridType == "G1") {
  TT <- TT.G1
} else if (details$gridType == "24-2") {
  TT <- TT.24.2
} else if (details$gridType == "practice") {
  TT <- tt.practice
} else if (details$gridType == "P-Total") {
  TT <- TT.PTotal
} else if (details$gridType == "P-Central10") {
  TT <- TT.PCentral10
} else if (details$gridType == "P-Central26") {
  TT <- TT.PCentral26
} else if (details$gridType == "P-Edge") {
  TT <- TT.PEdge
} else if (details$gridType == "P-Peripheral") {
  TT <- TT.PPeri
}

if (details$eye == "left") {
  TT <- grid.flip(TT)
  TT.PPeri <- grid.flip(TT.PPeri)
  TT.PCentral10 <- grid.flip(TT.PCentral10)
}

if (details$gridType == "P-Total") {
  z1 <- Zest242(eye=details$eye, primaryStartValue=PSV, gridType="P-Central26",outlierValue=7,minInterStimInterval=0,tt=TT.PCentral26,moveProjector = TRUE,retest=details$retest)
  tkdestroy(tt)
  graphics.off()
  details$fovea = FALSE
  z2 <- Zest242(eye=details$eye, primaryStartValue=PSV, gridType="P-Peripheral",outlierValue=7,minInterStimInterval=0,tt=TT.PPeri,moveProjector = TRUE,retest=details$retest)
  z <- combine(z1,z2)
} else {
  z <- Zest242(eye=details$eye, primaryStartValue=PSV, gridType=details$gridType,outlierValue=7,minInterStimInterval=0,tt=TT,moveProjector = TRUE,retest=details$retest)
}
  np <- c(np,sum(unlist(z$np),na.rm=TRUE)) # avg np 282 for tt=30
#}
terminate <- Sys.time()
tkdestroy(tt)
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
  details$comments <- paste(details$comments,comments,sep=".")
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



