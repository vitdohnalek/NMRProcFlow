# =============================================================================
# Spectral alignment: peak detection, CluPA, PTW warping, least-squares
# alignment, and manual shift.
# =============================================================================

#------------------------------
# Peak detection for spectra
#------------------------------
# Input parameters
#   - X: spectral dataset in matrix format in which each row contains a single sample
#   - nDivRange: size of a single small segment after division of spectra, Default value: 64
#   - baselineThresh: removal of all the peaks with intensity lower than this threshold, Default value: 50000
# Output parameters
#   - peak lists of the spectra
detectSpecPeaks <- function (X, nDivRange, scales=seq(1, 16, 2), baselineThresh, SNR.Th=-1, ProgressFile=NULL)
{
  LOGFILE <- file.path(dirname(ProgressFile),"clupa.out")
  if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,"Ralign1D:  BEGIN detectSpecPeaks");
  nFea <- ncol(X)
  nSamp <- nrow(X)
  noiseEsp <- 0.005
  if (SNR.Th < 0) SNR.Th <- max(scales) * 0.05
  if( !is.null(ProgressFile) ) init_counter(ProgressFile, nSamp)
  pList <- foreach(i=1:nSamp) %dopar% {
     myPeakRes <- NULL
     mySpec <- X[i, ]
     for (k in 1:length(nDivRange)) {
        divR <- nDivRange[k]
        for (j in 1:(trunc(nFea/divR) - 3)) {
           startR <- (j - 1) * divR + 1
           if (startR >= nFea)  startR <- nFea
           endR <- (j + 3) * divR
           if (endR > nFea) endR <- nFea
           xRange <- mySpec[startR:endR]
           xMean <- mean(xRange)
           xMedian <- median(xRange)
           if ((xMean == xMedian) || abs(xMean - xMedian)/((xMean + xMedian) * 2) < noiseEsp) next
           peakInfo <- MassSpecWavelet::peakDetectionCWT(mySpec[startR:endR], scales = scales, SNR.Th = SNR.Th)
           majorPeakInfo <- peakInfo$majorPeakInfo
           if (length(majorPeakInfo$peakIndex) > 0) myPeakRes <- c(myPeakRes, majorPeakInfo$peakIndex + startR - 1)
        }
     }
     plst <- list(myPeakRes)
     plst[[1]] <- sort(unique(plst[[1]]))
     plst[[1]] <- plst[[1]][which(mySpec[ plst[[1]] ] > baselineThresh)]
     plst[[1]] <- sort(plst[[1]])
     if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
     if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,paste("Ralign1D:  Spectrum",i," Nb peaks =",length(plst[[1]])));
     plst
  }
  if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,"Ralign1D:  END detectSpecPeaks");
  pList <-  simplify2array(pList)
  return(pList)
}


#------------------------------
# CluPA alignment for multiple spectra
#------------------------------
# Input parameters
#   - X: spectral dataset in the matrix format in which each row contains a single sample
#   - peakList: peak lists of the spectra
#   - refInd: index of the reference spectrum
#   - maxShift:  maximum number of the points for a shift step
# Output parameters
#   - aligned spectra: same format as input
dohCluster <- function (X, peakList, refInd = 1, maxShift = 50, acceptLostPeak = TRUE, ProgressFile=NULL)
{
  refSpec = X[refInd, ]
  if( !is.null(ProgressFile) ) init_counter(ProgressFile, nrow(X))
  Y <- foreach(tarInd=1:nrow(X), .combine=rbind) %dopar% {
     if (tarInd == refInd) { refSpec; } else {
     tarSpec <- X[tarInd, ]
     myPeakList <- c(peakList[[refInd]], peakList[[tarInd]])
     myPeakLabel <- double(length(myPeakList))
     myPeakLabel[1:length(peakList[[refInd]])] <- 1
     res <- hClustAlign(refSpec, tarSpec, myPeakList, myPeakLabel, 1, length(tarSpec), maxShift = maxShift, acceptLostPeak = TRUE)
     if( !is.null(ProgressFile) ) inc_counter(ProgressFile, tarInd)
     res$tarSpec
  }}
  return(Y)
}

#------------------------------
# Reference spectrum determination
#------------------------------
FindRef <- function (peakList)
{
    opts <- list(chunkSize=2)
    disS <- foreach(refInd = 1:length(peakList), .combine = 'rbind', .options.nws=opts) %dopar% {
        V <- rep(NA, length(peakList))
        for (tarInd in 1:length(peakList)) if (refInd != tarInd) {
            V[tarInd] = 0
            for (i in 1:length(peakList[[tarInd]])) V[tarInd] = V[tarInd] + min(abs(peakList[[tarInd]][i] - peakList[[refInd]]))
        }
        V
    }
    sumDis = double(length(peakList))
    for (refInd in 1:length(peakList)) {
        disS[refInd, refInd] = 0
        sumDis[refInd] = sum(disS[refInd, ])
    }
    orderSumdis = order(sumDis)
    return(list(refInd = orderSumdis[1], orderSpec = orderSumdis))
}

#------------------------------
# Spectra alignment - see https://cran.r-project.org/web/packages/speaq/vignettes/speaq.pdf
#------------------------------
# Input parameters
#   - data: n x p datamatrix
#   - nDivRange: size of a single small segment after division of the whole spectrum, Default value: 64
#   - reference: number of the spectrum reference; if NULL, automatic detection, Default value: NULL
#   - baselineThresh: removal of all the peaks with intensity lower than this threshold, Default value: 50000
# Output parameters
#   - Y: n x p datamatrix
CluPA <- function(data, reference=reference, nDivRange, scales = seq(1, 16, 2), baselineThresh,  SNR.Th = -1, maxShift=50, ProgressFile=NULL)
{
  LOGFILE <- NULL
  if( !is.null(ProgressFile) ) LOGFILE <- file.path(dirname(ProgressFile),conf$LOGFILE)
  if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,"Ralign1D:  BEGIN CluPA", mode="at");

  ## Peak picking
  if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,paste("Ralign1D:  --- Peak detection : nDivRange =",nDivRange));
  startTime <- proc.time()
  peakList <- detectSpecPeaks(X=data, nDivRange=nDivRange, scales=scales, baselineThresh=baselineThresh, SNR.Th = SNR.Th, ProgressFile=ProgressFile)
  endTime <- proc.time()
  if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,paste("Ralign1D:  --- Peak detection time: ",(endTime[3]-startTime[3])," sec"));
  if( !is.null(ProgressFile) ) unlink(ProgressFile)

  ## Reference spectrum determination
  if (reference == 0) {
     if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,"Ralign1D:  --- Find the spectrum reference...");
     startTime <- proc.time()
     resFindRef<- FindRef(peakList)
     refInd <- resFindRef$refInd
     endTime <- proc.time()
     if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,paste("Ralign1D:  --- The reference is: ",refInd));
     if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,paste("Ralign1D:  --- Finding time: ",(endTime[3]-startTime[3])," sec"));
  } else  {
     refInd=reference
  }
  ## Spectra alignment to the reference
  if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,paste("Ralign1D:  --- Spectra alignment to the reference: maxShift =",maxShift));
  startTime <- proc.time()
  Y <- dohCluster(data, peakList=peakList, refInd=refInd, maxShift=maxShift, acceptLostPeak, ProgressFile=ProgressFile)
  endTime <- proc.time()
  if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,paste("Ralign1D:  --- Spectra alignment time: ",(endTime[3]-startTime[3])," sec"));
  if( !is.null(ProgressFile) ) Write.LOG(LOGFILE,"Ralign1D:  END CluPA");

  ## Output
  return(Y)
}

#------------------------------
# LS : Alignment of the selected PPM ranges
#------------------------------
RAlign1D <- function(specMat, zone, RELDECAL=0.05, idxSref=0, Selected=NULL, fapodize=FALSE, ProgressFile=NULL)
{
   # Alignment of the PPM range
   NBPASS <- 3
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, NBPASS)
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   apodize <- ifelse(fapodize,1,0)
   decal <- round((i2-i1)*RELDECAL)
   for( n in 1:NBPASS) {
       ret <- align_segment(specMat$int, segment_shifts( specMat$int, idxSref, decal, i1-1, i2-1, Selected-1), i1-1, i2-1, apodize, Selected-1)
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, n)
   }

   return(specMat)
}

#------------------------------
# CluPA : Alignment of the selected PPM ranges
#------------------------------
RCluPA1D <- function(specMat, zonenoise, zone, resolution=0.02, SNR=3, idxSref=0, Selected=NULL, ProgressFile=NULL, nuc="1H")
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )

   # Noise estimation
   if (is.na(zonenoise)) {
       PPM_NOISE_AREA <- default_noise_range(nuc)
   } else {
      PPM_NOISE_AREA <- c(min(zonenoise), max(zonenoise))
   }
   idx_Noise <- c( length(which(specMat$ppm>PPM_NOISE_AREA[2])),(which(specMat$ppm<=PPM_NOISE_AREA[1])[1]) )
   Vref <- spec_ref(specMat$int)
   ynoise <- C_noise_estimation(Vref,idx_Noise[1],idx_Noise[2])

   # Parameters
   baselineThresh <- SNR*mean( C_noise_estimate(specMat$int, idx_Noise[1],idx_Noise[2], 1) )
   nDivRange <- max( round(resolution/specMat$dppm,0), 64 )
   maxshift <- min( round(0.01/specMat$dppm), round(nDivRange/4) )

   # Subpart of spectra
   if( is.null(Selected)) M<-specMat$int[, c(i1:i2) ] else  M<-specMat$int[Selected, c(i1:i2) ];

    M.aligned <- CluPA(M, reference=idxSref, nDivRange, scales = seq(1, 8, 2),
                          baselineThresh,  SNR.Th = 0.1, maxShift=maxshift, ProgressFile=ProgressFile)

   if( is.null(Selected)) specMat$int[ ,c(i1:i2)] <- M.aligned else specMat$int[ Selected,c(i1:i2)] <- M.aligned

   return(specMat)
}

#------------------------------
# PTW : Alignment of the selected PPM ranges
#------------------------------
RWarp1D <- function(specMat, zone, idxSref=0, warpcrit=c("WCC","RMS"), Selected=NULL)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )

   if( is.null(Selected)) M<-specMat$int[, c(i1:i2) ] else  M<-specMat$int[Selected, c(i1:i2) ];
   if (is.null(Selected)) nspec <- specMat$nspec else nspec <- length(Selected);

   if (idxSref==0 || ( !is.null(Selected) && !(idxSref %in% Selected) )) {
      V      <- bestref(M, optim.crit=warpcrit)
      refid  <- V$best.ref
   } else {
     refid <- idxSref
   }

   ref    <- M[refid, ]
   ssampl <- M[ c(1:nspec)[-refid], ]
   out    <- ptw(ref, ssampl, warp.type = "individual", mode = "forward", optim.crit=warpcrit)
   #out    <- ptw(ref, ssampl, warp.type = "global", mode = "forward", init.coef = c(0, 1, 0), optim.crit=warpcrit)
   M      <- out$warped.sample
   M[is.na(M)] <- 0
   if( is.null(Selected)) specMat$int[ c(1:nspec)[-refid],c(i1:i2)] <- M else specMat$int[ Selected[-refid],c(i1:i2)] <- M

   return(specMat)
}

#------------------------------
# Shift of the selected PPM ranges
#------------------------------
RShift1D <- function(specMat, zone, RELDECAL=0, Selected=NULL)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   di <- round(RELDECAL / specMat$dppm,0);
   j1 <- i1 - di
   j2 <- i2 - di
   icte <- ifelse(di<0, i1, i2)

   if( is.null(Selected) ) {
        M <- specMat$int[, c(i1:i2) ]
        for (k in 1:nrow(specMat$int)) specMat$int[ k, c(i1:i2) ] <- specMat$int[ k, icte ]
        specMat$int[, c(j1:j2) ] <- M
   } else  {
        M <- specMat$int[Selected, c(i1:i2) ]
        for (k in Selected) specMat$int[ k, c(i1:i2) ] <- specMat$int[ k, icte ]
        specMat$int[Selected, c(j1:j2) ] <- M
   }

   return(specMat)
}
