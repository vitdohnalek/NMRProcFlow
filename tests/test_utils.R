context("utils.R — helper functions")

# ---- default_noise_range ----

test_that("default_noise_range returns correct range for 1H", {
    expect_equal(default_noise_range("1H"), c(10.2, 10.5))
})

test_that("default_noise_range returns correct range for 13C", {
    expect_equal(default_noise_range("13C"), c(-20, -5))
})

test_that("default_noise_range falls back for unknown nucleus", {
    expect_equal(default_noise_range("31P"), c(10.2, 10.5))
})

test_that("default_noise_range falls back for empty string", {
    expect_equal(default_noise_range(""), c(10.2, 10.5))
})

# ---- trim ----

test_that("trim removes leading and trailing whitespace", {
    expect_equal(trim("  hello  "), "hello")
    expect_equal(trim("\thello\t"), "hello")
})

test_that("trim handles empty string", {
    expect_equal(trim(""), "")
})

test_that("trim leaves clean strings unchanged", {
    expect_equal(trim("hello"), "hello")
})

# ---- .N, .C, .toM ----

test_that(".N converts to numeric vector", {
    expect_equal(.N(c("1", "2.5", "3")), c(1, 2.5, 3))
})

test_that(".C converts factor to character vector", {
    f <- factor(c("a", "b", "c"))
    result <- .C(f)
    expect_equal(result, c("a", "b", "c"))
    expect_true(is.character(result))
})

test_that(".toM converts a matrix to numeric matrix preserving names", {
    # .toM calls as.numeric(x) which works on matrices but not data frames.
    # In the codebase, .toM is called on matrix objects, so test that path.
    M_in <- matrix(c("1", "2", "3", "4"), nrow = 2, ncol = 2)
    colnames(M_in) <- c("a", "b")
    rownames(M_in) <- c("r1", "r2")
    M <- .toM(M_in)
    expect_true(is.matrix(M))
    expect_true(is.numeric(M))
    expect_equal(M[1, 1], 1)
    expect_equal(M[2, 2], 4)
    expect_equal(colnames(M), c("a", "b"))
    expect_equal(rownames(M), c("r1", "r2"))
})

# ---- INI parsing round-trip ----

test_that("Write.INI and Parse.INI round-trip", {
    tmp <- tempfile(fileext = ".ini")
    on.exit(unlink(tmp))

    metalist <- list(PROCPARAMS = list(alpha = "10", beta = "TRUE", gamma = "hello"))
    Write.INI(tmp, metalist)

    result <- Parse.INI(tmp, section = "PROCPARAMS")
    expect_equal(result$alpha, 10)
    expect_equal(result$beta, TRUE)
    expect_equal(result$gamma, "hello")
})

# ---- Write.LOG ----

test_that("Write.LOG writes and appends", {
    tmp <- tempfile(fileext = ".log")
    on.exit(unlink(tmp))

    Write.LOG(tmp, "line one", mode = "wt")
    Write.LOG(tmp, "line two", mode = "at")

    lines <- readLines(tmp)
    expect_equal(lines, c("line one", "line two"))
})

# ---- counter functions ----

test_that("init/inc/get counter tracks progress", {
    tmpdir <- tempdir()
    pf <- file.path(tmpdir, "test_progress.log")
    on.exit({
        unlink(pf)
        unlink(file.path(tmpdir, ".out"), recursive = TRUE)
    })

    init_counter(pf, 5)
    result <- get_counter(pf)
    expect_equal(result$size, 5)
    expect_equal(result$value, 0)

    inc_counter(pf, 1)
    inc_counter(pf, 2)
    result <- get_counter(pf)
    expect_equal(result$value, 2)
})
