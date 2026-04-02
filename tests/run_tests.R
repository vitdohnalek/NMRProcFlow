#!/usr/bin/env Rscript
# =============================================================================
# Test runner for NMRProcFlow pure-R unit tests.
# Uses a lightweight base-R test framework (no testthat dependency).
#
# Usage:
#   Rscript tests/run_tests.R            # run all tests
#   Rscript tests/run_tests.R <file>     # run a single test file
#
# Dependencies: R + foreach + doParallel
# =============================================================================

# Resolve the directory of this script
.run_tests_dir <- (function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    }
    return("tests")
})()

# Source the helper (sets up test framework + loads modules)
source(file.path(.run_tests_dir, "helper_setup.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 1) {
    test_file <- args[1]
    if (!grepl("\\.R$", test_file)) test_file <- paste0(test_file, ".R")
    if (!file.exists(test_file)) {
        test_file <- file.path(.run_tests_dir, basename(test_file))
    }
    source(test_file)
} else {
    test_files <- sort(list.files(.run_tests_dir, pattern = "^test_.*\\.R$", full.names = TRUE))
    for (f in test_files) {
        source(f)
    }
}
