# =============================================================================
# RnmrTools.R — main loader
#
# This file sources the individual module files that contain all processing
# functions. Each module is self-contained and can be edited independently.
#
# Module layout (exec/lib/):
#   utils.R         — string helpers, type coercion, counters, INI/LOG I/O,
#                     default_noise_range()
#   stack.R         — undo history (push/pop/clean STACK)
#   metadata.R      — sample metadata generation from raw spectra archives
#   baseline.R      — baseline correction (Whittaker, airPLS, global, local, qNMR)
#   alignment.R     — spectral alignment (CluPA, PTW, least-squares, shift)
#   normalization.R — calibration and normalization (CSN, PQN)
#   processing.R    — denoising, zeroing, smoothing
#   bucketing.R     — spectral binning and macro-command validation
#   macro.R         — macro-command file replay (RProcCMD1D)
#   export.R        — bucket/SNR dataset generation and spectra data export
# =============================================================================

suppressMessages(library(Matrix))
suppressMessages(library(MASS))
suppressMessages(library(signal))
suppressMessages(library(ptw))
suppressMessages(library(speaq))
suppressMessages(library(rjson))

options(show.error.locations = TRUE)

# Resolve the directory containing this file so module paths work regardless
# of the caller's working directory.
# sys.frame(1)$ofile is set when source() is used; fall back to "exec" for
# cases where the file is evaluated in a different context (e.g. Rscript).
.RnmrTools_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) "exec"
)

source(file.path(.RnmrTools_dir, "lib", "utils.R"))
source(file.path(.RnmrTools_dir, "lib", "stack.R"))
source(file.path(.RnmrTools_dir, "lib", "metadata.R"))
source(file.path(.RnmrTools_dir, "lib", "baseline.R"))
source(file.path(.RnmrTools_dir, "lib", "alignment.R"))
source(file.path(.RnmrTools_dir, "lib", "normalization.R"))
source(file.path(.RnmrTools_dir, "lib", "processing.R"))
source(file.path(.RnmrTools_dir, "lib", "bucketing.R"))
source(file.path(.RnmrTools_dir, "lib", "macro.R"))
source(file.path(.RnmrTools_dir, "lib", "export.R"))
