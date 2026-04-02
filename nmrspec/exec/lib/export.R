# =============================================================================
# Export functions: bucket file parsing, bucket/SNR dataset generation,
# and raw spectra data export.
# =============================================================================

#----
# Generates the buckets table
#----
# Parse a bucket_list.in file into a tidy data frame.
# Columns returned: center, width, min, max, name
# Zero-width buckets are dropped; abs(width) guards against negative widths.
read_bucket_file <- function(bucketfile) {
    buckets <- read.table(bucketfile, header=F, sep="\t", stringsAsFactors=FALSE)
    buckets <- buckets[ buckets[,2] > 0, ]
    colnames(buckets) <- c("center", "width")
    buckets$min  <- buckets$center - 0.5 * abs(buckets$width)
    buckets$max  <- buckets$center + 0.5 * abs(buckets$width)
    buckets$name <- gsub("^(-?\\d+)", "B\\1",
                         gsub("\\.", "_",
                              gsub(" ", "", sprintf("%7.4f", buckets$center))))
    buckets
}

#----
get_Buckets_table <- function(bucketfile)
{
   outtable <- NULL
   if ( file.exists(bucketfile) ) {
      buckets <- read_bucket_file(bucketfile)
      outtable <- buckets[, c("name", "center", "min", "max", "width") ]
   }
   return(outtable)
}

#----
# Generates the buckets data set
#----
get_Buckets_dataset <- function(specMat, bucketfile, norm_method='CSN', zoneref=NA, YMAX=FALSE)
{
   outdata <- NULL
   if ( file.exists(bucketfile) ) {
      # Read the buckets
      buckets <- read_bucket_file(bucketfile)

      # get index of buckets' ranges
      buckets_m <- t(simplify2array(lapply( c( 1:(dim(buckets)[1]) ),  function(x){
             c(length(which(specMat$ppm>buckets[x,]$max)), length(which(specMat$ppm>buckets[x,]$min)))} )))

      # Takes the maximum intensity in the interval of each bucket rather than the integration
      if (YMAX) {
          buckets_IntVal <- NULL
          for (n in 1:specMat$nspec)
              buckets_IntVal <- rbind(buckets_IntVal, unlist(lapply(1:nrow(buckets_m), function(k){
                       max( specMat$int[n, buckets_m[k,1]:buckets_m[k,2] ] ) })))
      } else {
      # Integration
          buckets_IntVal <- C_all_buckets_integrate (specMat$int, buckets_m, 0)
          if (norm_method == 'CSN') {
              buckets_IntVal <- C_buckets_CSN_normalize( buckets_IntVal )
          }
          if (norm_method == 'PQN') {
              buckets_IntVal_CSN <- C_buckets_CSN_normalize( buckets_IntVal )
              bucVref_IntVal <- C_MedianSpec(buckets_IntVal_CSN)
              # For the quotient calculation, exclude columns where the reference is zero
              nonzero_cols <- bucVref_IntVal != 0
              if (!any(nonzero_cols)) {
                  warning("PQN: reference spectrum is entirely zero; skipping normalization")
              } else {
                  bucRatio <- sweep(buckets_IntVal_CSN[, nonzero_cols, drop=FALSE], 2,
                                    bucVref_IntVal[nonzero_cols], "/")
                  Coeff <- apply(bucRatio, 1, median)
                  bad_coeff <- Coeff == 0 | is.nan(Coeff) | is.infinite(Coeff)
                  if (any(bad_coeff)) {
                      warning(paste0("PQN: median quotient is zero/NaN/Inf for spectra: ",
                                     paste(which(bad_coeff), collapse=", "), "; leaving unnormalized"))
                      Coeff[bad_coeff] <- 1.0
                  }
                  # Apply Coeff to the full-width CSN matrix (not the filtered subset)
                  buckets_IntVal <- sweep(buckets_IntVal_CSN, 1, Coeff, "/")
              }
          }
      }

      # if supplied, integrate of all spectra within the PPM range of the reference signal
      if (! is.na(zoneref)) {
          istart <- length(which(specMat$ppm>max(zoneref)))
          iend <- length(which(specMat$ppm>min(zoneref)))
          Vref <- C_spectra_integrate (specMat$int, istart, iend)
          bad_ref <- Vref == 0 | is.nan(Vref) | is.infinite(Vref)
          if (any(bad_ref)) {
              warning(paste0("Reference signal normalization: integral is zero/NaN/Inf for spectra: ",
                             paste(which(bad_ref), collapse=", "), "; leaving unnormalized"))
              Vref[bad_ref] <- 1.0
          }
          buckets_IntVal <- buckets_IntVal/Vref
      }
      # Bucket names
      bucnames <- buckets$name

      # read samples
      samplesFile <- file.path(dirname(bucketfile),'samples.csv')
      samples <- read.table(samplesFile, header=F, sep=";", stringsAsFactors=FALSE)

      # read factors
      factorsFile <- file.path(dirname(bucketfile),"factors")
      factors <- read.table(factorsFile, header=F, sep=";", stringsAsFactors=FALSE)

      # Get Vendor and Pulse information
      Vendor <- toupper(readLines(file.path(dirname(bucketfile),'origin.txt'))[1])
      paramsFile <- file.path(dirname(bucketfile),'list_pars.csv')
      paramsDF <- read.table(paramsFile, header=T, sep=";", stringsAsFactors=FALSE)
      PULSE <- toupper((as.list(paramsDF[1,]))$PULSE)

      # Get the P15 parameter if Vendor == bruker and PULSE <=> cp (See Rnmr1D)
      is_cp <- length(grep(conf$CPREGEX, PULSE))>0
      if (Vendor == "BRUKER" && is_cp) {
          P15 <- paramsDF$P15
          outdata <- cbind( samples[, -1], P15, buckets_IntVal )
          colnames(outdata) <- c( factors[,2], 'P15', bucnames )
      } else {
          outdata <- cbind( samples[, -1], buckets_IntVal )
          colnames(outdata) <- c( factors[,2], bucnames )
      }
   }

   # return the data table
   return(outdata)
}

#----
# Generates the SNR dataset
#----
get_SNR_dataset <- function(specMat, bucketfile, zone_noise, ratio=TRUE)
{
   outdata <- NULL
   if ( file.exists(bucketfile) ) {
      # read samples
      samplesFile <- file.path(dirname(bucketfile),'samples.csv')
      samples <- read.table(samplesFile, header=F, sep=";", stringsAsFactors=FALSE)
      # read factors
      factorsFile <- file.path(dirname(bucketfile),"factors")
      factors <- read.table(factorsFile, header=F, sep=";", stringsAsFactors=FALSE)
      # Read the buckets
      buckets <- read_bucket_file(bucketfile)
      # get index of buckets' ranges
      buckets_m <- t(simplify2array(lapply( c( 1:(dim(buckets)[1]) ),
                     function(x) { c( length(which(specMat$ppm>buckets[x,]$max)), length(which(specMat$ppm>buckets[x,]$min)) ) }
                    )))
      # Compute Vnoise vector & Maxvals maxtrix
      i1 <- ifelse( max(zone_noise)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone_noise))) )
      i2 <- ifelse( min(zone_noise)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone_noise))[1] )
      flg <- 1
      Vnoise <- abs( C_noise_estimate(specMat$int, i1, i2, flg) )
      MaxVals <- C_maxval_buckets (specMat$int, buckets_m)
      # write the data table
      bucnames <- buckets$name
      if (ratio) {
         outdata <- cbind( samples[, -1], MaxVals/(2*Vnoise))
         colnames(outdata) <- c( factors[,2], bucnames )
      } else {
         outdata <- cbind( samples[, -1], Vnoise, MaxVals )
         colnames(outdata) <- c( factors[,2], 'Noise', bucnames )
      }

      #outdata <- data.frame(outdata, stringsAsFactors=FALSE)
   }
   return(outdata)
}


#----
# Generates the Spectra Data
#----
get_Spectra_Data <- function(specMat,samplesFile)
{
   # read samples
   samples <- read.table(samplesFile, header=F, sep=";", stringsAsFactors=FALSE)
   outdata <- cbind( specMat$ppm, t(specMat$int) )
   colnames(outdata) <- c( "ppm", samples[,1] )
   return(outdata)
}
