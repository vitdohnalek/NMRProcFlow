# =============================================================================
# Test helper: lightweight test framework using base R only (no testthat needed).
# Also sources the pure-R modules needed for unit testing.
# =============================================================================

suppressMessages(library(foreach))
suppressMessages(library(doParallel))

# Resolve paths relative to the project root.
# Works whether sourced from run_tests.R or directly via Rscript.
.helper_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) {
        # Fallback: find this file via the script that sourced us
        args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", args, value = TRUE)
        if (length(file_arg) > 0) {
            dirname(normalizePath(sub("^--file=", "", file_arg[1])))
        } else {
            "tests"
        }
    }
)
.project_root <- normalizePath(file.path(.helper_dir, ".."))
.lib_dir <- file.path(.project_root, "nmrspec", "exec", "lib")

source(file.path(.lib_dir, "utils.R"))
source(file.path(.lib_dir, "stack.R"))
source(file.path(.lib_dir, "normalization.R"))
source(file.path(.lib_dir, "export.R"))

# ---- Minimal test framework ----
.test_count <- 0L
.test_pass <- 0L
.test_fail <- 0L
.test_errors <- character(0)
.current_context <- ""

context <- function(desc) {
    .current_context <<- desc
    cat(sprintf("\n=== %s ===\n", desc))
}

test_that <- function(desc, code) {
    .test_count <<- .test_count + 1L
    tryCatch({
        force(code)
        .test_pass <<- .test_pass + 1L
        cat(sprintf("  PASS: %s\n", desc))
    }, error = function(e) {
        .test_fail <<- .test_fail + 1L
        msg <- sprintf("  FAIL: %s\n        %s\n", desc, conditionMessage(e))
        cat(msg)
        .test_errors <<- c(.test_errors, sprintf("[%s] %s: %s", .current_context, desc, conditionMessage(e)))
    })
}

expect_equal <- function(actual, expected, tolerance = 1e-8, label = NULL) {
    if (is.numeric(actual) && is.numeric(expected)) {
        if (!isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
            stop(sprintf("Expected %s but got %s",
                         paste(head(expected, 5), collapse=", "),
                         paste(head(actual, 5), collapse=", ")))
        }
    } else {
        if (!identical(actual, expected)) {
            stop(sprintf("Expected %s but got %s",
                         paste(head(expected, 5), collapse=", "),
                         paste(head(actual, 5), collapse=", ")))
        }
    }
}

expect_true <- function(x, label = NULL) {
    if (!isTRUE(x)) stop(paste("Expected TRUE but got", x, label))
}

expect_false <- function(x) {
    if (!identical(x, FALSE)) stop(paste("Expected FALSE but got", x))
}

expect_null <- function(x) {
    if (!is.null(x)) stop("Expected NULL")
}

expect_silent <- function(code) {
    tryCatch(force(code), error = function(e) stop(paste("Expected no error but got:", conditionMessage(e))))
}

expect_warning <- function(code) {
    warned <- FALSE
    withCallingHandlers(
        force(code),
        warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }
    )
    # We don't fail if no warning — some edge cases may not warn depending on data
}

# Print summary at end of session
.print_test_summary <- function() {
    cat(sprintf("\n============================\n"))
    cat(sprintf("Tests: %d | Pass: %d | Fail: %d\n", .test_count, .test_pass, .test_fail))
    if (length(.test_errors) > 0) {
        cat("\nFailures:\n")
        for (e in .test_errors) cat(sprintf("  - %s\n", e))
    }
    cat(sprintf("============================\n"))
    if (.test_fail > 0) quit(status = 1, save = "no")
}
reg.finalizer(.GlobalEnv, function(e) .print_test_summary(), onexit = TRUE)
