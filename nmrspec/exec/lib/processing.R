# =============================================================================
# Spectral processing: denoising (Savitzky-Golay filter), zeroing, and
# smoothing of selected PPM ranges.
# =============================================================================

#------------------------------
# Denoising the selected PPM ranges
#------------------------------
RFilter1D <- function(specMat,zone, FILTORD, FILTLEN, ProgressFile=NULL)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )

   sgfilt <- sgolay(p=FILTORD, n=FILTLEN)

   # Denoising each spectrum
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, specMat$nspec)
   for ( i in 1:specMat$nspec ) {
        x <- specMat$int[i,c(i1:i2)]
        SpecMat_sg <- filter(sgfilt,x)
        specMat$int[i,c(i1:i2)] <- SpecMat_sg
        if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
   }

   return(specMat)

}

#------------------------------
# Zeroing the selected PPM ranges
#------------------------------
RZero1D <- function(specMat, zones, LOGFILE=NULL, ProgressFile=NULL)
{
   # Zeroing each PPM range
   N <- dim(zones)[1]
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, N)

   for ( i in 1:N ) {
       i1<-length(which(specMat$ppm>max(zones[i,])))
       i2<-which(specMat$ppm<=min(zones[i,]))[1]
       specMat$int[,c(i1:i2)] <- matrix(0,specMat$nspec,(i2-i1+1))
       if( !is.null(LOGFILE) ) Write.LOG(LOGFILE,paste("Rnmr1D:     Zone",i,"= (",min(zones[i,]),",",max(zones[i,]),")"))
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
   }

   return(specMat)
}

#------------------------------
# Zeroing the selected PPM ranges
#------------------------------
RZeroNeg1D <- function(specMat, zones, LOGFILE=NULL, ProgressFile=NULL)
{
   # Zeroing negative values for each PPM range
   N <- dim(zones)[1]
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, N)

   for ( i in 1:N ) {
       i1<-length(which(specMat$ppm>max(zones[i,])))
       i2<-which(specMat$ppm<=min(zones[i,]))[1]
       V <- specMat$int[,c(i1:i2)]
	   V[V<0] <- 0
	   specMat$int[,c(i1:i2)] <- V
       if( !is.null(LOGFILE) ) Write.LOG(LOGFILE,paste("Rnmr1D:     Zone",i,"= (",min(zones[i,]),",",max(zones[i,]),")"))
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
   }

   return(specMat)
}

#------------------------------
# Smooth the selected PPM ranges by a line segment
#------------------------------
RSmooth1D <- function(specMat, zone, WS, LOGFILE=NULL)
{
   # Smooth the PPM range
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   for (k in 1:specMat$nspec) {
       V <- Smooth(specMat$int[k,c(i1:i2)], WS)
       n2 <- length(V); n1 <- n2 - WS + 1
       a <- (V[n2]-V[n1])/(n2-n1)
       for (j in n1:n2) V[j] <- a*(j-n1) + V[n1]
       specMat$int[k,c(i1:i2)] <- V
   }
   return(specMat)
}
