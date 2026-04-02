context("stack.R — undo history")

test_that("push_STACK creates versioned copies", {
    tmpdir <- tempfile("stack_test_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    writeLines("version0", file.path(tmpdir, "specs.pack"))
    push_STACK(tmpdir, "specs.pack", 0)

    expect_true(file.exists(file.path(tmpdir, "specs.pack.001")))
    expect_equal(readLines(file.path(tmpdir, "specs.pack.001")), "version0")
})

test_that("pop_STACK restores previous version", {
    tmpdir <- tempfile("stack_test_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    writeLines("version0", file.path(tmpdir, "specs.pack"))
    push_STACK(tmpdir, "specs.pack", 0)
    writeLines("version1", file.path(tmpdir, "specs.pack"))

    pop_STACK(tmpdir, "specs.pack", 1)
    expect_equal(readLines(file.path(tmpdir, "specs.pack")), "version0")
    expect_false(file.exists(file.path(tmpdir, "specs.pack.001")))
})

test_that("get_maxSTACKID returns correct max", {
    tmpdir <- tempfile("stack_test_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    writeLines("v0", file.path(tmpdir, "specs.pack"))
    expect_equal(get_maxSTACKID(tmpdir, "specs.pack"), 0)

    push_STACK(tmpdir, "specs.pack", 0)
    expect_equal(get_maxSTACKID(tmpdir, "specs.pack"), 1)

    push_STACK(tmpdir, "specs.pack", 1)
    expect_equal(get_maxSTACKID(tmpdir, "specs.pack"), 2)
})

test_that("clean_STACK removes all versioned files", {
    tmpdir <- tempfile("stack_test_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    writeLines("v0", file.path(tmpdir, "specs.pack"))
    push_STACK(tmpdir, "specs.pack", 0)
    push_STACK(tmpdir, "specs.pack", 1)
    push_STACK(tmpdir, "specs.pack", 2)

    clean_STACK(tmpdir, "specs.pack")

    expect_equal(get_maxSTACKID(tmpdir, "specs.pack"), 0)
    expect_true(file.exists(file.path(tmpdir, "specs.pack")))
})

test_that("push_STACK skips missing files without error", {
    tmpdir <- tempfile("stack_test_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    expect_silent(push_STACK(tmpdir, "nonexistent.file", 0))
})

test_that("pop_STACK skips missing versioned files without error", {
    tmpdir <- tempfile("stack_test_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    writeLines("current", file.path(tmpdir, "specs.pack"))
    expect_silent(pop_STACK(tmpdir, "specs.pack", 1))
    expect_equal(readLines(file.path(tmpdir, "specs.pack")), "current")
})
