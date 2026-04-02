# =============================================================================
# Baseline correction: Whittaker smoothing, airPLS, global/local/qNMR methods.
# =============================================================================

#------------------------------
#------------------------------
# airPLS
#------------------------------
WhittakerSmooth <- function(x,w,lambda,differences=1)
{
  x=matrix(x,nrow = 1, ncol=length(x))
  L=length(x)
  E=spMatrix(L,L,i=seq(1,L),j=seq(1,L),rep(1,L))
  D=methods::as(diff(E,1,differences),"CsparseMatrix") # 'dgCMatrix' move to 'CsparseMatrix'
  W=methods::as(spMatrix(L,L,i=seq(1,L),j=seq(1,L),w),"CsparseMatrix")
  background=solve((W+lambda*t(D)%*%D),t((w*x)));
  return(as.vector(background))
}

airPLS <- function(x,lambda=100, porder=1, itermax=8)
{
  x = as.vector(x)
  m = length(x)
  w = rep(1,m)
  i = 1
  repeat {
     z = WhittakerSmooth(x,w,lambda,porder)
     d = x-z
     sum_smaller = abs(sum(d[d<0]))
     if(sum_smaller<0.001*sum(abs(x))||i==itermax) break
     w[d>=0] = 0
     w[d<0] = exp(i*abs(d[d<0])/sum_smaller)
     w[1] = exp(i*max(d[d<0])/sum_smaller)
     w[m] = exp(i*max(d[d<0])/sum_smaller)
     i=i+1
  }
  return(z)
}

#------------------------------
# Global Baseline Correction
#------------------------------
RGbaseline1D <- function(specMat,PPM_NOISE_AREA, zone, WS, NEIGH, ProgressFile=NULL)
{
   NFAC <- 1.5
   NEGFAC <- 10
   NBPASS <- 2

   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   TD <- specMat$size

   # Baseline Estimation for each spectrum
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, specMat$nspec)
   BLList <- foreach(i=1:specMat$nspec, .combine=cbind) %dopar% {
       V <- specMat$int[i, ]
       if (WS>0) {
          sig <- fitdistr(specMat$int[i, length(which(specMat$ppm>PPM_NOISE_AREA[2])):(which(specMat$ppm<=PPM_NOISE_AREA[1])[1])], "normal")$estimate[2]
          mv <- simplify2array(lapply( c(3:61), function(x) { mean(abs(specMat$int[i, ((x-1)*TD/64):(x*TD/64)])); }))
          mmoy<-min(mv)
          if (min(specMat$int[i, ])< -NEGFAC*mmoy) {
              Vthreshold <- max(specMat$int[ i, specMat$int[i, ]/mmoy < -NEGFAC ])
              specMat$int[ i, specMat$int[i, ]/mmoy < -NEGFAC ] <- Vthreshold
          }
          if (NBPASS==1) {
              BL <- C_Estime_LB (V, i1, i2, WS, NEIGH, NFAC*sig)
          } else {
              BL <- 0*rep(1:length(V))
              for( n in 1:NBPASS) {
                 # Estimation of Baseline
                 BLn <- C_Estime_LB (V, i1, i2, WS, NEIGH, NFAC*sig)
                 V <- V - BLn
                 BL <- BL + BLn
              }
          }
       } else {
          BL <- ptw::asysm(V, lambda = 1e+10, p = 0.05, eps = 1e-8, maxit = 42)
       }
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
       BL
   }

   # Baseline Correction for each spectrum
   for ( i in 1:specMat$nspec ) {
       mv <- simplify2array(lapply( c(3:61), function(x) { mean(abs(specMat$int[i, ((x-1)*TD/64):(x*TD/64)])); }))
       mmoy<-min(mv)
       if (min(specMat$int[i, ])< -NEGFAC*mmoy) {
           Vthreshold <- max(specMat$int[ i, specMat$int[i, ]/mmoy < -NEGFAC ])
           specMat$int[ i, specMat$int[i, ]/mmoy < -NEGFAC ] <- Vthreshold
       }
       if( is.null(dim(BLList)) ) { BL <- BLList; } else { BL <- BLList[,i]; }
       specMat$int[i,c(i1:i2)] <- specMat$int[i,c(i1:i2)] - BL[c(i1:i2)]
   }

   return(specMat)
}

#------------------------------
# Local Baseline Correction (deprecated)
#------------------------------
RBaseline1D <- function(specMat,PPM_NOISE_AREA, zone, WINDOWSIZE, ProgressFile=NULL)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )

   SMOOTHSIZE <- round(WINDOWSIZE/2)
   ALPHA <- 0.2

   # Ajust the window size parameter
   n <- i2-i1+1
   ws <- WINDOWSIZE
   dws <- 0
   signdws <- ifelse ( n/ws>round(n/ws), 1, -1 )
   dmin <- 1
   if (n>WINDOWSIZE) {
      repeat {
          d <- abs(n/(ws+signdws*dws)-round(n/(ws+signdws*dws)))
          if ( d>dmin ) { dws <- dws - signdws; break }
          dmin <- d; dws <- dws + signdws;
      }
      WINDOWSIZE <- ws+signdws*dws
      n2 <- round(n/(ws+signdws*dws))*(ws+signdws*dws)
      if (n2<n) i2 <- i2 - (n-n2)
   } else {
      WINDOWSIZE <- n
   }

   # Baseline Estimation for each spectrum
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, specMat$nspec)
   BLList <- foreach(i=1:specMat$nspec, .combine=cbind) %dopar% {
       specSig <- fitdistr(specMat$int[i,length(which(specMat$ppm>PPM_NOISE_AREA[2])):(which(specMat$ppm<=PPM_NOISE_AREA[1])[1])], "normal")$estimate[2]
       x <- specMat$int[i,c(i1:i2)]
       xmat <- matrix(x, nrow=WINDOWSIZE)
       ymin <- apply(xmat, 2, min) + 1.2*specSig
       y1 <- rep(ymin, each = WINDOWSIZE)
       y2 <- Smooth(y1,SMOOTHSIZE)
       BL <- simplify2array(lapply(c(1:length(x)), function(k) { min(y1[k], y2[k]); }))
       BL <- lowpass1(BL, ALPHA)
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
       BL
   }

   # Baseline Correction for each spectrum
   for ( i in 1:specMat$nspec ) {
       if( is.null(dim(BLList)) ) { BL <- BLList; } else { BL <- BLList[,i]; }
       specMat$int[i,c(i1:i2)] <- specMat$int[i,c(i1:i2)] - BL
   }

   return(specMat)
}

#------------------------------
# q-NMR Baseline Correction
#------------------------------
Rqnmrbc1D <- function(specMat, PPM_NOISE_AREA, zone, ProgressFile=NULL)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   n <- i2-i1+1
   NLOOP <- 5
   dN <- round(0.00075/specMat$dppm)
   CSIG <- 5

   # Baseline Estimation for each spectrum
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, specMat$nspec)
   inoise <- length(which(specMat$ppm>PPM_NOISE_AREA[2])):(which(specMat$ppm<=PPM_NOISE_AREA[1])[1])
   BLList <- foreach(i=1:specMat$nspec, .combine=cbind) %dopar% {
       specSig <- fitdistr(specMat$int[i,inoise], "normal")$estimate[2]
       x <- specMat$int[i,c(i1:i2)]
       bc <- 0*rep(1:n)
       for (l in 1:NLOOP) {
           bci <- C_GlobSeg(x, dN, CSIG*specSig)
           x <- x - bci
           bc <- bc + bci
       }
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
       bc
  }

   # Baseline Correction for each spectrum
   for ( i in 1:specMat$nspec ) {
       if( is.null(dim(BLList)) ) { BL <- BLList; } else { BL <- BLList[,i]; }
       V <- specMat$int[i,c(i1:i2)] - BL
       specMat$int[i,c(i1:i2)] <- V
   }

   return(specMat)
}

#------------------------------
# airPLS : Local Baseline Correction
#------------------------------
RairPLSbc1D <- function(specMat, zone, clambda, porder=1, ProgressFile=NULL)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   n <- i2-i1+1
   cmax <- switch(porder, 6, 7, 8)

   #lambda <- ifelse (clambda==cmax, 5, 10^(cmax-clambda) )
   lambda <- ifelse( clambda>0, cmax-clambda, abs(clambda) )
   lambda <- 10^max(lambda,1)

   # Baseline Estimation for each spectrum
   if( !is.null(ProgressFile) ) init_counter(ProgressFile, specMat$nspec)
   BLList <- foreach(i=1:specMat$nspec, .combine=cbind) %dopar% {
       bc <- airPLS(specMat$int[i,c(i1:i2)], lambda, porder=porder)
       if( !is.null(ProgressFile) ) inc_counter(ProgressFile, i)
       bc
   }

   # Baseline Correction for each spectrum
   for ( i in 1:specMat$nspec ) {
       if( is.null(dim(BLList)) ) { BL <- BLList; } else { BL <- BLList[,i]; }
       V <- specMat$int[i,c(i1:i2)] - BL
       specMat$int[i,c(i1:i2)] <- V
   }

   return(specMat)
}
