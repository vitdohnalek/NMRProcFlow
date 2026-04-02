# =============================================================================
# Macro-command processing: replay a sequence of spectral processing commands
# from a macro-command file.
# =============================================================================

# ------------------------------------
# Process the Macro-commands file
# ------------------------------------
RProcCMD1D <- function(specMat, specParamsDF, CMDTEXT, NCPU=1, LOGFILE=NULL, ProgressFile=NULL)
{
   lbALIGN <- 'align'
   lbGBASELINE <- 'gbaseline'
   lbBASELINE <- 'baseline'
   lbQNMRBL <- 'qnmrbline'
   lbAIRPLS <- 'airpls'
   lbBIN <- 'binning'
   lbCALIB <- 'calibration'
   lbNORM <- 'normalisation'
   lbFILTER <- 'denoising'
   lbWARP <- 'warp'
   lbCLUPA <- 'clupa'
   lbSHIFT <- 'shift'
   lbBUCKET <- 'bucket'
   lbZERO <- 'zero'
   lbZERONEG <- 'zeroneg'
   lbSMOOTH <- 'smooth'
   EOL <- 'EOL'

   SI <- (as.list(specParamsDF[1,]))$SI

   CMDTEXT <- CMDTEXT[ grep( "^[^ ]", CMDTEXT ) ]
   CMD <- CMDTEXT[ grep( "^[^#]", CMDTEXT ) ]
   CMD <- gsub("^ ", "", gsub(" $", "", gsub(" +", ";", CMD)))

   samples <- read.table( 'samples.csv', header=F, sep=";", stringsAsFactors=FALSE)

   specMat$fWriteSpec <- FALSE

   while ( length(CMD)>0 && CMD[1] != EOL ) {

      cmdLine <- CMD[1]
      cmdPars <- unlist(strsplit(cmdLine[1],";"))
      cmdName <- cmdPars[1]

      repeat {
          if (cmdName == lbCALIB) {
              params <- as.numeric(cmdPars[c(-1,-7)])
              if (length(params)>=3) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 PPMREF <- params[3]
                 PPM_NOISE <- ifelse( length(params)>4, c( min(params[4:5]), max(params[4:5]) ), c( 10.2, 10.5 ) )
                 CALIBTYPE <- ifelse(  cmdPars[7] %in% c('s','d'), cmdPars[7], 's' )
                 Write.LOG(LOGFILE, paste0("Rnmr1D:  Calibration: PPM REF =",PPMREF,", Zone Ref = (",PPMRANGE[1],",",PPMRANGE[2],"), Type = ",CALIBTYPE));
                 registerDoParallel(cores=NCPU)
                 specMat <- RCalib1D(specMat, PPM_NOISE, PPMRANGE, PPMREF, CALIBTYPE, ProgressFile=ProgressFile)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              if (length(params)==1) {
                 ppmshift <- params[1]
                 Write.LOG(LOGFILE, paste0("Rnmr1D:  Calibration: PPM shift =",ppmshift));
                 specMat$ppm_min <- specMat$ppm_min + ppmshift
                 specMat$ppm_max <- specMat$ppm_max + ppmshift
                 specMat$ppm <- rev(seq(from=specMat$ppm_min, to=specMat$ppm_max, by=specMat$dppm))
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbNORM) {
              params <- cmdPars[-1]
              if (length(params)==2) {
                 # Legacy 2-parameter format: "normalisation <ppm_min> <ppm_max>"
                 # Method is not stored in this format; default to CSN and warn.
                 params <- as.numeric(params)
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Normalisation (legacy format, defaulting to CSN): Zone Ref = (",PPMRANGE[1],",",PPMRANGE[2],")"))
                 warning("Macro replay: legacy 2-parameter normalisation format does not encode the method; replaying as CSN. Re-save the macro to preserve the original method.")
                 registerDoParallel(cores=NCPU)
                 specMat <- RNorm1D(specMat, norm_method='CSN', zones=matrix(PPMRANGE,nrow=1, ncol=2))
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              if (length(params)==1) {
                 SPECTRAL_NORM_METH <- params[1]
                 CMD <- CMD[-1]
                 zones <- NULL
                 while(CMD[1] != EOL) {
                    zones <- rbind(zones, as.numeric(unlist(strsplit(CMD[1],";"))))
                    CMD <- CMD[-1]
                 }
                 Write.LOG(LOGFILE,"Rnmr1D:  Normalisation of the Intensities based on the selected PPM ranges...")
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Method =",SPECTRAL_NORM_METH))
                 registerDoParallel(cores=NCPU)
                 specMat <- RNorm1D(specMat, norm_method=SPECTRAL_NORM_METH, zones=zones)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbGBASELINE || cmdName == lbBASELINE) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)==6) {
                 PPM_NOISE <- c( min(params[1:2]), max(params[1:2]) )
                 PPMRANGE <- c( min(params[3:4]), max(params[3:4]) )
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Baseline Correction: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"));
                 registerDoParallel(cores=NCPU)
                 if (cmdName == lbBASELINE) {
                     if (params[5]>0) { # to be compliant with version<1.1.4
                        BCMETH <- params[5]
                        WSFAC <- (7-params[6])/4
                        WINDOWSIZE <- round(WSFAC*SI/(BCMETH*64))
                     } else {
                        WINDOWSIZE <- round(( 1/2^(params[6]-2) )*(SI/64))
                     }
                     Write.LOG(LOGFILE,paste0("Rnmr1D:     Type=Local - Window Size = ",WINDOWSIZE));
                     specMat <- RBaseline1D(specMat,PPM_NOISE, PPMRANGE, WINDOWSIZE, ProgressFile=ProgressFile)
                 } else {
                     WS <- params[5]
                     NEIGH <- params[6]
                     Write.LOG(LOGFILE,paste0("Rnmr1D:     Type=Global - Smoothing Parameter=",WS," - Window Size=",NEIGH));
                     specMat <- RGbaseline1D(specMat,PPM_NOISE, PPMRANGE, WS, NEIGH, ProgressFile=ProgressFile)
                 }
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbQNMRBL) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)==4) {
                 PPM_NOISE <- c( min(params[1:2]), max(params[1:2]) )
                 PPMRANGE <- c( min(params[3:4]), max(params[3:4]) )
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Baseline Correction: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Type=q-NMR"))
                 registerDoParallel(cores=NCPU)
                 specMat <- Rqnmrbc1D(specMat,PPM_NOISE, PPMRANGE, ProgressFile=ProgressFile)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbAIRPLS) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)>=3) {
                 porder <- 1
                 if (length(params)==4) porder <- params[4]
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 LAMBDA <- params[3]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Baseline Correction: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Type=airPLS, lambda=",LAMBDA, ", order=",porder))
                 registerDoParallel(cores=NCPU)
                 specMat <- RairPLSbc1D(specMat, PPMRANGE, LAMBDA, porder=porder, ProgressFile=ProgressFile)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }

          if (cmdName == lbFILTER) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)==4) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 FORDER <- params[3]
                 FLENGTH <- params[4]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Denoising: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"));
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Filter Order=",FORDER," - Filter Length=",FLENGTH));
                 registerDoParallel(cores=NCPU)
                 specMat <- RFilter1D(specMat,PPMRANGE, FORDER, FLENGTH, ProgressFile=ProgressFile)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbALIGN) {
              params <- as.numeric(cmdPars[-1])
              Selected <- NULL
              if (length(params)==6 && (params[5]<2 || params[6])) {
                 level <- unique(samples[ order(as.character(samples[, params[5]+1])), params[5]+1 ])[params[6]]
                 Selected <- .N(rownames(samples[ samples[, params[5]+1]==level, ]))
              }
              if (length(params)>=4) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 RELDECAL= params[3]
                 idxSref=params[4]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Alignment: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Rel. Shift Max.=",RELDECAL," - Reference=",idxSref))
                 specMat <- RAlign1D(specMat, PPMRANGE, RELDECAL, idxSref, Selected=Selected, fapodize=FALSE, ProgressFile=ProgressFile)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbWARP) {
              params <- cmdPars[-1]
              Selected <- NULL
              if (length(params)==6 && (.N(params[5])<2 || .N(params[6]))) {
                 level <- unique(samples[ order(as.character(samples[, .N(params)[5]+1])), .N(params)[5]+1 ])[.N(params[6])]
                 Selected <- .N(rownames(samples[ samples[, .N(params[5])+1]==level, ]))
              }
              if (length(params)>=4) {
                 PPMRANGE <- c( min(.N(params[1:2])), max(.N(params[1:2])) )
                 idxSref=.N(params[3])
                 warpcrit=.C(params[4])
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Alignment: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Parametric Time Warping Method - Reference=",idxSref," - Optim. Crit=",warpcrit))
                 specMat <- RWarp1D(specMat, PPMRANGE, idxSref, warpcrit, Selected=Selected)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbCLUPA) {
              params <- as.numeric(cmdPars[-1])
              Selected <- NULL
              if (length(params)==9 && (params[8]<2 || params[9])) {
                 level <- unique(samples[ order(as.character(samples[, params[8]+1])), params[8]+1 ])[params[9]]
                 Selected <- .N(rownames(samples[ samples[, params[8]+1]==level, ]))
              }
              if (length(params)>=7) {
                 PPM_NOISE <- c( min(params[1:2]), max(params[1:2]) )
                 PPMRANGE <- c( min(params[3:4]), max(params[3:4]) )
                 RESOL= params[5]
                 SNR= params[6]
                 idxSref=params[7]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Alignment: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     CluPA - Resolution =",RESOL," - SNR threshold=",SNR, " - Reference=",idxSref))
                 registerDoParallel(cores=NCPU)
                 specMat <- RCluPA1D(specMat, PPM_NOISE, PPMRANGE, RESOL, SNR, idxSref, Selected=Selected, ProgressFile=ProgressFile)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbSHIFT) {
              params <- as.numeric(cmdPars[-1])
              Selected <- NULL
              if (length(params)==5 && (params[4]<2 || params[5])) {
                 level <- unique(samples[ order(as.character(samples[, params[4]+1])), params[4]+1 ])[params[5]]
                 Selected <- .N(rownames(samples[ samples[, params[4]+1]==level, ]))
              }
              if (length(params)>=3) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 RELDECAL= params[3]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Shift: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Shift value =",RELDECAL))
                 specMat <- RShift1D(specMat, PPMRANGE, RELDECAL, Selected=Selected)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbZERO) {
              CMD <- CMD[-1]
              zones2 <- NULL
              while(CMD[1] != EOL) {
                  zones2 <- rbind(zones2, as.numeric(unlist(strsplit(CMD[1],";"))))
                  CMD <- CMD[-1]
              }
              Write.LOG(LOGFILE,"Rnmr1D:  Zeroing the selected PPM ranges ...")
              specMat <- RZero1D(specMat, zones2, LOGFILE=LOGFILE, ProgressFile=ProgressFile)
              specMat$fWriteSpec <- TRUE
              CMD <- CMD[-1]
              break
          }
          if (cmdName == lbZERONEG) {
              CMD <- CMD[-1]
              zones2 <- NULL
              while(CMD[1] != EOL) {
                  zones2 <- rbind(zones2, as.numeric(unlist(strsplit(CMD[1],";"))))
                  CMD <- CMD[-1]
              }
              Write.LOG(LOGFILE,"Rnmr1D:  Zeroing negative values for the selected PPM ranges ...")
              specMat <- RZeroNeg1D(specMat, zones2, LOGFILE=LOGFILE, ProgressFile=ProgressFile)
              specMat$fWriteSpec <- TRUE
              CMD <- CMD[-1]
              break
          }
          if (cmdName == lbSMOOTH) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)>=3) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 WS <- params[3]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Smooth: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Window size =",WS))
                 specMat <- RSmooth1D(specMat, PPMRANGE, WS, LOGFILE=LOGFILE)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
                 break
              }
          }
          if (cmdName == lbBUCKET) {
              if ( !( length(cmdPars) >= 6 && cmdPars[2] %in% c('aibin','erva','unif') ) &&
                   !( length(cmdPars) == 2 && cmdPars[2] %in% c('vsb') ) ) {
                 CMD <- CMD[-1]
                 break;
              }
              Write.LOG(LOGFILE,"Rnmr1D: \nRnmr1D:  Bucketing the selected PPM ranges ...")
              CMD <- CMD[-1]
              zones <- NULL
              while(CMD[1] != EOL) {
                  zones <- rbind(zones, as.numeric(unlist(strsplit(CMD[1],";"))))
                  CMD <- CMD[-1]
              }
              if ( cmdPars[2] %in% c('aibin','erva','unif') ) {
                  params <- as.numeric(cmdPars[-c(1:2)])
                  PPM_NOISE <- c( min(params[1:2]), max(params[1:2]) )
                  resol <- params[3]; snr <- params[4];
                  Write.LOG(LOGFILE,paste0("Rnmr1D:     ",toupper(cmdPars[2])," - Resolution=",resol," - SNR threshold=",snr))
              } else {
                  PPM_NOISE <- NULL
                  resol <- 0; snr <- 0;
                  Write.LOG(LOGFILE,paste0("Rnmr1D:     ",toupper(cmdPars[2])))
              }
              registerDoParallel(cores=NCPU)
              RBucket1D(specMat, cmdPars[2], resol, snr, zones, PPM_NOISE, LOGFILE=LOGFILE, ProgressFile=ProgressFile)
              CMD <- CMD[-1]
              break
          }
          CMD <- CMD[-1]
          break
      }
      gc()
   }
   return(specMat)
}
