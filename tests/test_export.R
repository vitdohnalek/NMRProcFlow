context("export.R — bucket file parsing")

test_that("read_bucket_file parses a valid bucket file", {
    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp))

    writeLines(c(
        "1.5000\t0.1000",
        "2.5000\t0.2000",
        "3.5000\t0.0500"
    ), tmp)

    result <- read_bucket_file(tmp)

    expect_equal(nrow(result), 3)
    expect_equal(colnames(result), c("center", "width", "min", "max", "name"))
    expect_equal(result$center, c(1.5, 2.5, 3.5))
    expect_equal(result$width, c(0.1, 0.2, 0.05))
    expect_equal(result$min[1], 1.5 - 0.05)
    expect_equal(result$max[1], 1.5 + 0.05)
})

test_that("read_bucket_file drops zero-width buckets", {
    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp))

    writeLines(c(
        "1.5000\t0.1000",
        "2.5000\t0.0000",
        "3.5000\t0.0500"
    ), tmp)

    result <- read_bucket_file(tmp)
    expect_equal(nrow(result), 2)
    expect_equal(result$center, c(1.5, 3.5))
})

test_that("read_bucket_file generates valid bucket names", {
    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp))

    writeLines(c(
        "1.5000\t0.1000",
        "-0.3000\t0.0500"
    ), tmp)

    result <- read_bucket_file(tmp)

    expect_true(all(grepl("^B", result$name)))
    expect_false(any(grepl("\\.", result$name)))
    expect_true(grepl("^B-", result$name[2]))
})

test_that("read_bucket_file drops negative-width buckets (filtered by > 0)", {
    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp))

    writeLines(c(
        "5.0000\t0.2000",
        "6.0000\t-0.1000"
    ), tmp)

    result <- read_bucket_file(tmp)
    # Negative width row is dropped by the > 0 filter
    expect_equal(nrow(result), 1)
    expect_equal(result$center[1], 5.0)
})

test_that("get_Buckets_table returns NULL for missing file", {
    expect_null(get_Buckets_table("/nonexistent/path/buckets.txt"))
})

test_that("get_Buckets_table returns correct columns", {
    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp))

    writeLines(c(
        "1.5000\t0.1000",
        "2.5000\t0.2000"
    ), tmp)

    result <- get_Buckets_table(tmp)
    expect_equal(colnames(result), c("name", "center", "min", "max", "width"))
    expect_equal(nrow(result), 2)
})
