# =============================================================================
# Stack management: undo history for successive versions of the binary spectra
# matrix (specs.pack file).
# =============================================================================

### Management of the historic of processings, (i.e. the sucessive versions of the matrix of the binary spectra, i.e. the specs.pack file)
#   push : Save the input (add) on the top of the stack (i.e. the current matrix of the binary spectra, the logfile)
#   pop  : Remove the last input from the stack (i.e. the last saved matrix of the binary spectra on the top of the list of files)  and replace by the previous one

get_maxSTACKID <- function(DATADIR, FILENAME)
{
   lstFiles <- list.files(path=DATADIR, full.names=FALSE, recursive = FALSE)
   if (length( grep(paste0(FILENAME,'.[0-9]'), lstFiles) )==0) return(0)
   return(max(as.numeric(gsub(paste0(FILENAME,'.'), '', lstFiles[grep(paste0(FILENAME,'.[0-9]'), lstFiles)]))))
}

push_STACK <- function (DATADIR, listfiles, STACKID)
{
   for (f in 1:length(listfiles)) {
      F0 <- file.path(DATADIR,listfiles[f])
      if (! file.exists(F0)) next
      FN <- file.path(DATADIR,paste0(listfiles[f],'.',sprintf("%03d",STACKID+1)))
      file.copy(F0, FN, overwrite = TRUE)
   }
}

pop_STACK <- function (DATADIR, listfiles, STACKID)
{
   for (f in 1:length(listfiles)) {
      FN <- file.path(DATADIR,paste0(listfiles[f],'.',sprintf("%03d",STACKID)))
      F0 <- file.path(DATADIR,listfiles[f])
      if (! file.exists(FN)) next
      file.copy(FN, F0, overwrite = TRUE)
      unlink(FN)
   }
}

clean_STACK <- function(DATADIR, listfiles)
{
   for (f in 1:length(listfiles))
      file.remove( file.path(DATADIR, dir(path=DATADIR ,pattern=paste0(listfiles[f],'.0*'))) )
}
