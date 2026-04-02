# =============================================================================
# Normalization and calibration of spectra.
# =============================================================================

#------------------------------
# Calibration ot the PPM Scale
#------------------------------
RCalib1D <- function(specMat, PPM_NOISE_AREA, zoneref, ppmref, type='s', ProgressFile=NULL)
{
   i1<-length(which(specMat$ppm>max(zoneref)))
   i2<-which(specMat$ppm<=min(zoneref))[1]
   PPM_MIN <- -1000
   PPM_MAX <- 1000
   N <- round((specMat$ppm_max-specMat$ppm_min)/specMat$dppm)
   DMIN <- round(10*N/65535)

   # Compute the shift of each spectrum
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, specMat$nspec)
   Tdecal <- foreach::foreach(i=1:specMat$nspec, .combine=c) %dopar% {
       if (type=='d') {
           V <- order(specMat$int[i, i1:i2], decreasing=T)
           k <- 2; while(abs(V[k]-V[k-1])<DMIN) k <- k+1
           i0 <- i1 + round(0.5*(V[1]+V[k])) - 1
       } else {
           i0 <- i1 + which(specMat$int[i, i1:i2]==max(specMat$int[i, i1:i2])) - 1
       }
       ppm0 <- specMat$ppm_max - (i0-1)*specMat$dppm
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
       return(ppm0 - ppmref)
   }

   # Compute the new full ppm range
   PPM_MIN <- max(PPM_MIN, specMat$ppm_min-Tdecal)
   PPM_MAX <- min(PPM_MAX, specMat$ppm_max-Tdecal)

   # PPM calibration of each spectrum
   N <- length(seq(from=PPM_MIN, to=PPM_MAX, by=specMat$dppm))
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, specMat$nspec)
   M <- foreach::foreach(i=1:specMat$nspec, .combine=rbind) %dopar% {
       ppm <- specMat$ppm - Tdecal[i]
       P <- ppm>PPM_MIN & ppm<=PPM_MAX
       V <- specMat$int[i,P]
       if (length(V)<N) { V <- c(V, rep(0,N-length(V)) ) }
       if (length(V)>N) { V <- V[1:N] }
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
       return(V)
   }
   specMat$int <- M
   specMat$ppm_min <- PPM_MIN
   specMat$ppm_max <- PPM_MAX
   specMat$ppm <- rev(seq(from=specMat$ppm_min, to=specMat$ppm_max, by=specMat$dppm))
   return(specMat)
}

#------------------------------
# Normalisation of the Intensities
#------------------------------
RNorm1D <- function(specMat, norm_method, zones)
{
   N <- dim(zones)[1]
   # Use sequential foreach for single zone to avoid requiring a registered backend
   if (N > 1) {
      registerDoParallel(cores=2)
      `%norm_op%` <- `%dopar%`
   } else {
      `%norm_op%` <- `%do%`
   }

   # Helper: integrate a single zone across all spectra (trapezoidal rule)
   .integrate_zone <- function(int_mat, ppm, zone) {
       i1 <- length(which(ppm > max(zone)))
       i2 <- which(ppm <= min(zone))[1]
       # Guard: zone outside spectral range or too narrow
       if (is.na(i1) || is.na(i2) || i1 == 0 || i2 <= i1) {
           warning(paste0("Normalization: zone (", min(zone), ", ", max(zone),
                          ") is outside the spectral range or too narrow; returning zeros"))
           return(rep(0, nrow(int_mat)))
       }
       simplify2array(lapply(1:nrow(int_mat), function(x) {
           if (i2 == i1 + 1) {
               # Only two boundary points, no interior — simple trapezoid
               0.5*(int_mat[x, i1] + int_mat[x, i2])
           } else {
               0.5*(int_mat[x, i1] + int_mat[x, i2]) + sum(int_mat[x, (i1+1):(i2-1)])
           }
       }))
   }

   # Helper: guard a coefficient vector against zero/NaN/Inf
   .guard_coeff <- function(coeff) {
       bad <- coeff == 0 | is.nan(coeff) | is.infinite(coeff)
       if (any(bad)) {
           warning(paste0("Normalization: coefficient is zero/NaN/Inf for spectra: ",
                          paste(which(bad), collapse=", "), "; leaving those spectra unnormalized"))
           coeff[bad] <- 1.0
       }
       coeff
   }

   if (norm_method=='CSN') {
      # Integrate each zone and sum across zones
      SUM <- foreach(i=1:N, .combine='+') %norm_op% {
          .integrate_zone(specMat$int, specMat$ppm, zones[i,])
      }
      meanSUM <- mean(SUM)
      if (meanSUM == 0 || is.nan(meanSUM)) {
          warning("CSN: mean normalization coefficient is zero/NaN; skipping normalization")
          return(specMat)
      }
      COEFF <- .guard_coeff(SUM / meanSUM)
   }

   if (norm_method=='PQN') {
      # Step 1: Initial CSN pre-normalization (required per Dieterle et al. 2006)
      SUM <- foreach(i=1:N, .combine='+') %norm_op% {
          .integrate_zone(specMat$int, specMat$ppm, zones[i,])
      }
      meanSUM <- mean(SUM)
      if (meanSUM == 0 || is.nan(meanSUM)) {
          warning("PQN: CSN pre-normalization mean is zero/NaN; skipping normalization")
          return(specMat)
      }
      CSN_COEFF <- .guard_coeff(SUM / meanSUM)
      CSN_INT <- specMat$int / CSN_COEFF

      # Steps 2-4: PQN on the CSN-normalized intensities
      SUBMAT <- foreach(i=1:N, .combine=cbind) %norm_op% {
          i1 <- length(which(specMat$ppm > max(zones[i,])))
          i2 <- which(specMat$ppm <= min(zones[i,]))[1]
          t(simplify2array(lapply(1:specMat$nspec, function(x) { CSN_INT[x, i1:i2] })))
      }
      V <- apply(SUBMAT, 2, median)
      SUBMAT <- SUBMAT[, V != 0, drop=FALSE]
      V <- V[V != 0]
      if (length(V) == 0) {
          warning("PQN: all reference column medians are zero; skipping normalization")
          return(specMat)
      }
      MQ <- t(t(SUBMAT) / V)
      PQN_COEFF <- apply(MQ, 1, median)
      COEFF <- .guard_coeff(CSN_COEFF * PQN_COEFF)
   }

   # Apply coefficient: divide each spectrum (row) by its coefficient
   specMat$int <- specMat$int / COEFF
   return(specMat)
}
