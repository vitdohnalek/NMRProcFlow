context("normalization.R — CSN and PQN")

# Helper: build a minimal specMat-like list for testing.
# ppm runs from ppm_max (index 1) down to ppm_min (last index), matching NMR convention.
# Use zones strictly inside (ppm_min, ppm_max) to ensure i1 >= 1 and i2 <= size.
make_specMat <- function(int_matrix, ppm_min = -0.5, ppm_max = 11) {
    nspec <- nrow(int_matrix)
    size <- ncol(int_matrix)
    dppm <- (ppm_max - ppm_min) / (size - 1)
    ppm <- rev(seq(from = ppm_min, to = ppm_max, by = dppm))
    if (length(ppm) > size) ppm <- ppm[1:size]
    if (length(ppm) < size) ppm <- c(ppm, rep(ppm[length(ppm)], size - length(ppm)))
    list(
        int = int_matrix,
        ppm = ppm,
        nspec = nspec,
        size = size,
        ppm_min = ppm_min,
        ppm_max = ppm_max,
        dppm = dppm
    )
}

# ============================================================================
# CSN (Constant Sum Normalization)
# ============================================================================

test_that("CSN normalizes spectra to equal total intensity", {
    M <- matrix(0, nrow = 3, ncol = 100)
    M[1, ] <- rep(10, 100)
    M[2, ] <- rep(20, 100)
    M[3, ] <- rep(30, 100)

    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    result <- RNorm1D(sm, norm_method = "CSN", zones = zones)

    # CSN equalizes zone integrals (trapezoidal), not raw rowSums.
    # For constant-value spectra, all values should converge to the mean (20).
    expect_equal(result$int[1, 50], result$int[2, 50], tolerance = 1e-6)
    expect_equal(result$int[2, 50], result$int[3, 50], tolerance = 1e-6)
})

test_that("CSN with identical spectra leaves them unchanged", {
    M <- matrix(5, nrow = 3, ncol = 50)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    result <- RNorm1D(sm, norm_method = "CSN", zones = zones)
    expect_equal(result$int, M)
})

test_that("CSN known-answer: spec1=10, spec2=30 converge to same value", {
    M <- matrix(c(rep(10, 50), rep(30, 50)), nrow = 2, byrow = TRUE)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    result <- RNorm1D(sm, norm_method = "CSN", zones = zones)

    # After CSN, both constant spectra should have the same value everywhere.
    # COEFF[1] = 10/20 = 0.5, COEFF[2] = 30/20 = 1.5
    # Result[1] = 10/0.5 = 20, Result[2] = 30/1.5 = 20
    # Check interior points (boundary points may differ slightly due to trapezoidal half-weight)
    expect_true(all(abs(result$int[1, 5:45] - result$int[2, 5:45]) < 1e-6))
    expect_true(abs(result$int[1, 25] - 20) < 1e-6)
})

test_that("CSN with single spectrum leaves it unchanged", {
    set.seed(42)
    M <- matrix(runif(100, 1, 10), nrow = 1)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    result <- RNorm1D(sm, norm_method = "CSN", zones = zones)
    expect_equal(result$int, M)
})

test_that("CSN guards against all-zero spectrum (no Inf/NaN)", {
    M <- matrix(c(rep(10, 50), rep(0, 50)), nrow = 2, byrow = TRUE)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    expect_warning(result <- RNorm1D(sm, norm_method = "CSN", zones = zones))
    expect_false(any(is.nan(result$int)))
    expect_false(any(is.infinite(result$int)))
})

test_that("CSN with entirely zero spectra skips normalization", {
    M <- matrix(0, nrow = 3, ncol = 50)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    expect_warning(result <- RNorm1D(sm, norm_method = "CSN", zones = zones))
    expect_equal(result$int, M)
})

test_that("CSN with multiple zones sums across zones", {
    M <- matrix(c(rep(10, 100), rep(20, 100)), nrow = 2, byrow = TRUE)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 5, 5, 10), nrow = 2, ncol = 2, byrow = TRUE)

    result <- RNorm1D(sm, norm_method = "CSN", zones = zones)

    # Both spectra should converge to the same values
    expect_equal(result$int[1, 50], result$int[2, 50], tolerance = 1e-6)
})

# ============================================================================
# PQN (Probabilistic Quotient Normalization)
# ============================================================================

test_that("PQN performs CSN pre-step then equalizes scaled spectra (NORM-03 regression)", {
    M <- matrix(0, nrow = 4, ncol = 100)
    M[1, ] <- seq(1, 100)
    M[2, ] <- seq(2, 200, by = 2)
    M[3, ] <- seq(3, 300, by = 3)
    M[4, ] <- seq(0.5, 50, by = 0.5)

    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    result <- RNorm1D(sm, norm_method = "PQN", zones = zones)

    # Spectra with identical shape (just scaled) should converge
    ratio_12 <- result$int[1, ] / result$int[2, ]
    expect_true(sd(ratio_12) / mean(ratio_12) < 0.01)
})

test_that("PQN with identical spectra leaves them unchanged", {
    M <- matrix(7, nrow = 3, ncol = 50)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    result <- RNorm1D(sm, norm_method = "PQN", zones = zones)
    expect_equal(result$int, M)
})

test_that("PQN with single spectrum leaves it unchanged", {
    set.seed(42)
    M <- matrix(runif(100, 1, 10), nrow = 1)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    result <- RNorm1D(sm, norm_method = "PQN", zones = zones)
    expect_equal(result$int, M)
})

test_that("PQN with all-zero spectra skips normalization", {
    M <- matrix(0, nrow = 3, ncol = 50)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    expect_warning(result <- RNorm1D(sm, norm_method = "PQN", zones = zones))
    expect_equal(result$int, M)
})

test_that("PQN guards against zero/NaN/Inf coefficients", {
    M <- matrix(c(rep(10, 50), rep(0, 50), rep(5, 50)), nrow = 3, byrow = TRUE)
    sm <- make_specMat(M)
    zones <- matrix(c(0, 10), nrow = 1, ncol = 2)

    expect_warning(result <- RNorm1D(sm, norm_method = "PQN", zones = zones))
    expect_false(any(is.nan(result$int)))
    expect_false(any(is.infinite(result$int)))
    # Zero spectrum should remain zero
    expect_true(all(result$int[2, ] == 0))
})

# ============================================================================
# Edge cases: zone boundaries
# ============================================================================

test_that("normalization with zone outside spectral range warns", {
    M <- matrix(10, nrow = 2, ncol = 50)
    sm <- make_specMat(M)
    zones <- matrix(c(20, 30), nrow = 1, ncol = 2)

    expect_warning(result <- RNorm1D(sm, norm_method = "CSN", zones = zones))
})

test_that("normalization with narrow zone does not produce NaN/Inf", {
    M <- matrix(c(rep(10, 50), rep(20, 50)), nrow = 2, byrow = TRUE)
    sm <- make_specMat(M)
    zones <- matrix(c(4.9, 5.1), nrow = 1, ncol = 2)

    result <- RNorm1D(sm, norm_method = "CSN", zones = zones)
    expect_false(any(is.nan(result$int)))
    expect_false(any(is.infinite(result$int)))
})
