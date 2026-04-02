# =============================================================================
# Metadata generation: unzip raw data, parse sample files, write metadata.
# =============================================================================

# -----
# Generate the 'samples.csv' & 'factors' files from the list of raw spectra
# -----

# RAWDIR : the directory containing the uploaded files (unzipped ZIP and samples files)
#   samples.txt : the uploaded Sample file
# DATADIR : where data will be stored and used by NMRViewer
#   samples.csv : the Sample file used by NMRViewer
#   factors : file given the list of the severals factors used by NMRViewer
generate_Metadata_File <- function(RawZip, DATADIR, procParams)
{
   # Unzip the archive
   RAWDIR <- dirname(RawZip)
   ext <- tolower(gsub("^.*\\.", "", RawZip))
   if (ext=='7z') {
       # Use system2 with explicit argument vector to prevent shell injection via filename
       system2("7zr", args=c("x", "-y", paste0("-o", RAWDIR), RawZip))
   } else {
       unzip(RawZip, files = NULL, list = FALSE, overwrite = TRUE,  junkpaths = FALSE, exdir = RAWDIR, unzip = "internal",   setTimes = FALSE)
   }

   # if a file  of samples was uploaded along with the ZIP file
   samples <- NULL
   SampleFile <- file.path(RAWDIR,'samples.txt')
   if (file.exists(SampleFile)) {
        min_col <- 2
        if ( procParams$VENDOR=="bruker" || procParams$VENDOR=="rs2d" ) min_col <- 4
        samples <- read.table(SampleFile, sep="\t", header=T,stringsAsFactors=FALSE)
        if (ncol(samples)< min_col) samples <- NULL
   }
   metadata <- generateMetadata(RAWDIR, procParams, samples)

   if (length(metadata$ERRORLIST)>0) {
      write.table(metadata$ERRORLIST, file=file.path(DATADIR,'errorlist.csv'), sep=';', row.names=F, col.names=F, quote=F)
   }
   OKRAW <- 0

   if ( (class(metadata$rawids)=="matrix" && dim(metadata$rawids)[1]>1) || length(metadata$rawids)>0 ) {
      write.table(metadata$samples, file=file.path(DATADIR,'samples.csv'), sep=';', row.names=F, col.names=F, quote=F)
      write.table(metadata$rawids, file=file.path(DATADIR,'rawids.csv'), sep=';', row.names=F, col.names=F, quote=F)
      write.table(metadata$factors, file=file.path(DATADIR,'factors'), sep=';', row.names=F, col.names=F, quote=F)
      OKRAW <- 1
   }

   return(OKRAW)
}
