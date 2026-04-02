# =============================================================================
# Utility functions: string helpers, type coercion, progress counters,
# INI file I/O, and log writing.
# =============================================================================

# Default noise PPM range for noise estimation, chosen per nucleus.
# The region must be signal-free in typical metabolomics spectra.
#   1H : 10.2-10.5 ppm  — beyond aromatic/aldehyde signals
#   13C : -20 to -5 ppm — sub-zero region, no metabolite 13C signals exist here
#         (do NOT use 200-220 ppm: that overlaps with carbonyls and may exceed PPM_MAX_13C=200)
default_noise_range <- function(nuc) {
    switch(nuc,
        "1H"  = c(10.2, 10.5),
        "13C" = c(-20, -5),
        c(10.2, 10.5)        # fallback for unrecognised nuclei
    )
}

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
.N <- function(x) { as.numeric(as.vector(x)) }
.C <- function(x) { as.vector(x) }
.toM <- function(x) { M <- matrix(as.numeric(x), nrow=dim(x)[1], ncol=dim(x)[2]);
                     colnames(M) <- colnames(x); rownames(M) <- rownames(x); M }

# Init counter
init_counter <- function(ProgressFile, n) {
    fh<-file(ProgressFile,"wt"); writeLines(as.character(n), fh); close(fh)
    dirout <- file.path(dirname(ProgressFile),".out")
    if( file.exists(dirout) ) unlink(dirout, recursive=TRUE)
    dir.create(dirout)
}

# Put i to 1
inc_counter <- function(ProgressFile, i) {
   fileout <- file.path(dirname(ProgressFile),".out",paste0(i,".out"))
   if (! file.exists(fileout)) {
       fl<-file(fileout,"wt"); writeLines("1", fl); close(fl)
   }
}

# Read the counter
get_counter <- function(ProgressFile) {
    n <- as.numeric(readLines(ProgressFile))
    c <- length(list.files(path=file.path(dirname(ProgressFile),".out"), pattern="*.out", all.files=FALSE, full.names=FALSE))
    return( list(value=c, size=n))
}

### Write within the 'INI.file' file  with the INI format, the 'metalist' list, a list  of lists
#   EXCLU being a list of keys to exclude to the writting
Write.INI <- function(INI.file, metalist, EXCLU=c())
{
   INI.list <- c()
   for (i in 1:length(metalist)) {
      section <- names(metalist)[i]
      INI.list <- c( INI.list, paste0('[',section,']') )
      M<-unlist(metalist[[section]])
      for ( k in 1:length(M) ) {
          if ( names(M[k]) %in% EXCLU ) next
          INI.list <- c( INI.list, paste0(names(M[k]),'=',M[k],sep="") )
      }
      INI.list <- c( INI.list, '' )
   }
   write.table(INI.list, file=INI.file, sep='', row.names=F, col.names=F, quote=F)
}

### Parse the section 'section' within the 'INI.file' file
#   Get the INI.list as an initial list to add or replace the couple of values (key=value)
Parse.INI <- function(INI.file, INI.list=list(), section="PROCPARAMS")
{
   connection <- file(INI.file)
   Lines  <- readLines(connection)
   close(connection)

   Lines <- chartr("[]", "==", Lines)  # change section headers

   connection <- textConnection(Lines)
   d <- read.table(connection, as.is = TRUE, sep = "=", fill = TRUE)
   close(connection)

   L <- d$V1 == ""                    # location of section breaks
   d <- subset(transform(d, V3 = V2[which(L)[cumsum(L)]])[1:3], V1 != "")
   d <- d[d$V3 == section,]

   #INI.list <- list()
   for( i in 1:dim(d)[1] ) {
        key <- d$V1[i]
        val <- d$V2[i]
        if (! is.na(suppressWarnings(as.numeric(val)))) {
            INI.list[[key]] <- as.numeric(val)
            next
        }
        if (! is.na(suppressWarnings(as.logical(val)))) {
            INI.list[[key]] <- as.logical(val)
            next
        }
        INI.list[[key]] <- val
   }
   return(INI.list)
}

### Write within the 'LOG.file' file, the 'textline' text
#   mode : can be either 'at' for 'appending' mode or 'wt' for 'writing' mode
Write.LOG <- function(LOG.file, textline="", mode="at")
{
   fileLog<-file(LOG.file,mode)
   writeLines(textline, fileLog)
   close(fileLog)
}
