# =============================================================================
# Bucketing: spectral binning algorithms (intelligent, uniform, ERVA, VSB)
# and macro-command file validation/processing.
# =============================================================================

#------------------------------
# Bucket : Apply the bucketing based on the 'Algo' algorithm with the resolution 'resol'.
# Then elinate buckets with a SNR under the threshold given by 'snr'
#------------------------------
RBucket1D <- function(specMat, Algo, resol, snr, zones, zonenoise, LOGFILE=NULL, ProgressFile=NULL)
{
   BUCKET_LIST  <- 'bucket_list.in'
   BUC.filename <- 'SpecBuckets.txt'
   BUC.cmd <- 'SpecBucCmd.lst'
   NUC <- readLines('nuc.txt')
   wrtCMD <- FALSE
   idx <- 1:specMat$nspec

   # Limit size of buckets
   MAXBUCKETS<-2000
   NOISE_FAC <- 3

   if (Algo %in% c('aibin','unif','erva')) {
      # Noise estimation
      if (is.na(zonenoise)) {
          PPM_NOISE_AREA <- default_noise_range(NUC)
      } else {
         PPM_NOISE_AREA <- c(min(zonenoise), max(zonenoise))
      }
      idx_Noise <- c( length(which(specMat$ppm>PPM_NOISE_AREA[2])),(which(specMat$ppm<=PPM_NOISE_AREA[1])[1]) )
      Vref <- spec_ref(specMat$int)
      ynoise <- C_noise_estimation(Vref,idx_Noise[1],idx_Noise[2])
      Vnoise <- abs( C_noise_estimate(specMat$int, idx_Noise[1],idx_Noise[2], 1) )
      if (wrtCMD) {
          Write.LOG(BUC.cmd,sprintf("#\n# Bucketing - Method: %s, Noise Zone = (%f,%f), Resolution = %f, SNR =  %d\n#",
               toupper(algo), PPM_NOISE_AREA[1], PPM_NOISE_AREA[2], resol, snr))
          Write.LOG(BUC.cmd,sprintf("bucket %s %f %f %f %d",Algo, PPM_NOISE_AREA[1], PPM_NOISE_AREA[2], resol, snr))
      }
   }

   if (Algo %in% c('aibin')) {
      bdata <- list()
      bdata$ynoise <- ynoise
      bdata$vnoise <- NULL
      bdata$inoise_start <- idx_Noise[1]
      bdata$inoise_end <- idx_Noise[2]
      bdata$R <- resol
      bdata$dppm <- specMat$dppm
      bdata$noise_fac <- NOISE_FAC
      bdata$VREF <- 1
      if (NUC=="13C") {
         bdata$noise_fac <- 2
         bdata$bin_fac <- 0.1
         bdata$peaknoise_rate <- 5
         bdata$BUCMIN <- 0.05
      } else {
         bdata$noise_fac <- NOISE_FAC
         bdata$bin_fac <- 0.5
         bdata$peaknoise_rate <- 15
         bdata$BUCMIN <- 0.003
      }
   }

   if (Algo %in% c('erva')) {
      bdata <- list()
      bdata$bucketsize <- resol
      bdata$noise_fac <- 1
      bdata$dppm <- specMat$dppm
      bdata$ppm_min <- specMat$ppm_min
      bdata$BUCMIN <- 0.001
      # if CP sequence then withdrraw the first spectra from the kinetics (tc<TCmin)
      if (!is.null(specParamsDF) && 'P15' %in% colnames(specParamsDF)) {
          TCmin <- ifelse(!is.null(procParams) && !is.null(procParams$TCmin), procParams$TCmin, 100)
          idx <- which(specParamsDF$P15>=TCmin)
          if( !is.null(LOGFILE) ) Write.LOG(LOGFILE,paste("Rnmr1D:     CP sequence: TCmin =",TCmin))
      }
   }

   if (Algo=='vsb' && wrtCMD) {
      Write.LOG(BUC.cmd,sprintf("bucket %s",Algo), mode="at")
   }

   # For each PPM range
   buckets_zones <- NULL
   N <- dim(zones)[1]
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, N)
   buckets_zones <- foreach(i=1:N, .combine=rbind) %dopar% {
       if (wrtCMD) Write.LOG(BUC.cmd,paste(min(zones[i,]), max(zones[i,])))
       i2<-which(specMat$ppm<=min(zones[i,]))[1]
       i1<-length(which(specMat$ppm>max(zones[i,])))
       if (Algo=='aibin') {
          Mbuc <- matrix(, nrow = MAXBUCKETS, ncol = 2)
          Mbuc[] <- 0
          buckets_m <- C_aibin_buckets(specMat$int, Mbuc, Vref, bdata, i1, i2)
       }
       if (Algo=='erva') {
          Mbuc <- matrix(, nrow = MAXBUCKETS, ncol = 2)
          Mbuc[] <- 0
          buckets_m <- C_erva_buckets(specMat$int[idx, ], Mbuc, Vref, bdata, i1, i2)
          V <- abs(specMat$ppm[buckets_m[,2]] - specMat$ppm[buckets_m[,1]])
          buckets_m <- buckets_m[ which(V>5*specMat$dppm), ]
       }
       if (Algo=='unif') {
          seq_buc <- seq(i1, i2, round(resol/specMat$dppm))
          n_bucs <- length(seq_buc) - 1
          buckets_m <- cbind ( seq_buc[1:n_bucs], seq_buc[2:(n_bucs+1)])
       }
       if (Algo=='vsb') {
          buckets_m <- matrix( c( i1, i2 ), nrow=1, ncol=2, byrow=T )
       }
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
       # Keep only the buckets for which the SNR average is greater than 'snr'
       if (nrow(buckets_m)>1) {
          MaxVals <- C_maxval_buckets (specMat$int, buckets_m)
          bucsel <- which( apply(t(MaxVals/(2*Vnoise)),1,stats::quantile)[4,]>snr )
          buckets_m <- buckets_m[ bucsel, ]
       }
       if( !is.null(LOGFILE) )
          Write.LOG(LOGFILE,paste("Rnmr1D:     Zone",i,"= (",min(zones[i,]),",",max(zones[i,]),"), Nb Buckets =",nrow(buckets_m)))
       cbind( specMat$ppm[buckets_m[,1]], specMat$ppm[buckets_m[,2]] )
   }
   if (wrtCMD) Write.LOG(BUC.cmd,"EOL\n", mode="at")

   if( !is.null(LOGFILE) ) Write.LOG(LOGFILE,paste("Rnmr1D:     Total Buckets =",nrow(buckets_zones)))

   if( buckets_zones[1,1]>buckets_zones[1,2] )  {  colnames(buckets_zones) <- c('max','min') }
                                          else  {  colnames(buckets_zones) <- c('min','max') }
   f_append <- ifelse (file.exists(BUC.filename), TRUE, FALSE )

   # The bucket zones files
   write.table(buckets_zones, file=BUC.filename, append=f_append, sep="\t", row.names=F, col.names=!f_append, quote=F)
   buclist <- cbind( 0.5*(buckets_zones[,1]+buckets_zones[,2]), abs(buckets_zones[,2]-buckets_zones[,1]) )
   write.table(buclist, file=BUCKET_LIST, append=f_append, sep="\t", row.names=F, col.names=F, quote=F)
   gc()
}

#------------------------------
# Check if the macro-command file (CMD.filename) is compliant with the allowed commands
#------------------------------
check_MacroCmdFile <- function(CMD.filename) {
   ret <- 1
   allowKW <- c( 'align', 'warp', 'clupa', 'shift', 'gbaseline', 'baseline', 'qnmrbline', 'airpls', 'binning', 'calibration', 'normalisation', 'denoising', 'bucket', 'zero', 'zeroneg', 'smooth', 'EOL' )

   tryCatch({
      # Read the macrocommand file
      CMDTEXT <- gsub("\t", "", readLines(CMD.filename))
      CMDTEXT <- CMDTEXT[ grep( "^[^ ]", CMDTEXT ) ]
      CMD <- CMDTEXT[ grep( "^[^#]", CMDTEXT ) ]
      CMD <- gsub("^ ", "", gsub(" $", "", gsub(" +", ";", CMD)))
      L <- unique(sort(gsub(";.*$","", CMD)))
      L <- L[ grep( "^[^0-9-]", L)]
      ret <- ifelse( sum(L %in% allowKW)==length(L), 1, 0 )
   }, error=function(e) {
       ret <- 0
   })
   return(ret)
}
