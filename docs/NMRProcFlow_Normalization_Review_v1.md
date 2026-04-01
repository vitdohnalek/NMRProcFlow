# NMRProcFlow Code Review — Normalization Deep-Dive & General Issues

**Repository:** NMRProcFlow-master (v1.2.12)
**Review date:** 2026-03-19
**Scope:** Full codebase analysis with primary focus on normalization correctness
**Key files examined:**
- `nmrspec/exec/RnmrTools.R` — core processing algorithms
- `nmrspec/exec/libspec/libCspec.cpp` — C++ performance-critical routines
- `nmrspec/exec/Rcorr1D` — preprocessing orchestration script
- `nmrspec/Rsrc/Proc4.R` — Shiny server logic for processing
- `nmrspec/Rsrc/ui_procparams.R` — UI parameter definitions
- `nmrspec/Rsrc/utils.R` — data export and workbook generation
- `nmrspec/exec/Rbuc1D` — bucketing module

---

## Table of Contents

1. [Normalization Architecture Overview](#1-normalization-architecture-overview)
2. [NORM-01: Two Divergent Normalization Pipelines](#norm-01-two-divergent-normalization-pipelines)
3. [NORM-02: CSN Algorithm Mismatch Between Pipelines](#norm-02-csn-algorithm-mismatch-between-pipelines)
4. [NORM-03: PQN Missing Required CSN Pre-Normalization (Path A)](#norm-03-pqn-missing-required-csn-pre-normalization-path-a)
5. [NORM-04: Biased Median in C_MedianSpec](#norm-04-biased-median-in-c_medianspec)
6. [NORM-05: Division-by-Zero — No Guards Anywhere](#norm-05-division-by-zero--no-guards-anywhere)
7. [NORM-06: Silent Triple-Layer Normalization Stacking](#norm-06-silent-triple-layer-normalization-stacking)
8. [NORM-07: Parallelization Bug in RNorm1D](#norm-07-parallelization-bug-in-rnorm1d)
9. [NORM-08: Hardcoded Macro Replay Forces CSN](#norm-08-hardcoded-macro-replay-forces-csn)
10. [NORM-09: Inconsistent Parameter Naming](#norm-09-inconsistent-parameter-naming)
11. [General Issues](#general-issues)
12. [Proposed Fix Priority](#proposed-fix-priority)

---

## 1. Normalization Architecture Overview

Before diving into individual bugs, it is essential to understand that normalization in NMRProcFlow is **not a single operation** — it occurs in two completely separate code paths, implemented independently, using different algorithms for the same named methods.

### Path A — Spectral-Level Normalization (Preprocessing)

```
User selects: Processing tab → "Normalisation" → CSN or PQN
                                ↓
Proc4.R:get_CorrCmd()  →  sets procParams$NORM_METH
                                ↓
Rcorr1D (line ~160)    →  reads zones from zones1_list.in
                                ↓
RNorm1D()              →  modifies specMat$int directly
  (RnmrTools.R:421)       writes updated specs.pack to disk
```

This operates on the **full spectral intensity matrix** and permanently alters the stored spectra.

### Path B — Bucket-Level Normalization (Data Export)

```
User selects: Data Export tab → Normalization Method → None/CSN/PQN
                                ↓
utils.R:get_Data_matrix()  →  calls get_Buckets_dataset()
                                ↓
get_Buckets_dataset()      →  integrates buckets first
  (RnmrTools.R:1329)          then applies C_buckets_CSN_normalize (C++)
                               or PQN using C_MedianSpec (C++)
```

This operates on **already-integrated bucket values** and only affects the exported data, not the stored spectra.

### Path C — Macro Replay Normalization

```
Macro command file → Rnmr1D_DoProc() → parses "normalisation" keyword
  (RnmrTools.R:1050)                    calls RNorm1D() (same as Path A)
```

This is a third entry point that replays saved processing commands.

**These three paths are the source of most normalization problems documented below.**

---

## NORM-01: Two Divergent Normalization Pipelines

| | |
|---|---|
| **Severity** | Critical |
| **Files** | `RnmrTools.R:421–458` (Path A), `RnmrTools.R:1329–1402` + `libCspec.cpp:1075–1092` (Path B) |
| **Impact** | Silent double-normalization; different results from same named method |

### Problem

Nothing prevents a user from applying normalization during preprocessing (Path A via the Processing tab) and then selecting normalization again at export (Path B via the Data Export tab). The UI has no state tracking and no warning.

When this happens, the data is **double-normalized**. Worse, even if the user selects the same method both times (e.g., CSN + CSN), the results are **not equivalent to applying it once**, because the two CSN implementations use different algorithms (see NORM-02).

### How a user triggers this

1. Upload spectra → Processing tab → select "Normalisation" → CSN → Process
2. Go to Data Export tab → Normalization Method dropdown shows "None" (default) but user changes to "CSN" → Export

The exported data is now CSN-on-CSN with two different CSN algorithms.

### Recommendation

- Track whether spectral-level normalization has been applied (flag in the session state).
- When a user selects export normalization and spectral normalization was already applied, either disable the export option or display a prominent warning: *"Normalization was already applied during preprocessing. Applying it again will double-normalize your data."*
- Consider removing one of the two normalization points entirely and consolidating into a single, well-tested path.

---

## NORM-02: CSN Algorithm Mismatch Between Pipelines

| | |
|---|---|
| **Severity** | Critical |
| **Files** | `RnmrTools.R:426–436` vs. `libCspec.cpp:1075–1092` |
| **Impact** | "CSN" produces fundamentally different results depending on where it's applied |

### Path A — CSN in `RNorm1D` (R code)

```r
# RnmrTools.R, lines 428-435
SUM <- foreach(i=1:N, .combine='+') %dopar% {
    i1 <- length(which(specMat$ppm > max(zones[i,])))
    i2 <- which(specMat$ppm <= min(zones[i,]))[1]
    simplify2array(lapply(1:specMat$nspec, function(x) {
        0.5*(specMat$int[x, i1] + specMat$int[x, i2]) +
            sum(specMat$int[x, (i1+1):(i2-1)])
    }))
}
COEFF <- SUM / mean(SUM)           # <-- relative to mean
```

**What it does:** Computes the trapezoidal integral for each spectrum over the selected zones, then divides each integral by the mean of all integrals. The coefficient is a ratio centered on the group mean. A spectrum with exactly average intensity gets `COEFF = 1.0` and is unchanged.

**Mathematical form:** `COEFF[k] = SUM[k] / mean(SUM)`, then `spectrum[k] = spectrum[k] / COEFF[k]`

### Path B — CSN in `C_buckets_CSN_normalize` (C++ code)

```cpp
// libCspec.cpp, lines 1085-1089
for (k = 0; k < n_specs; k++) {
    sumS = 0.0;
    for (m = 0; m < n_bucs; m++) sumS += buckets(k, m);
    for (m = 0; m < n_bucs; m++) M(k, m) = 100000.0 * buckets(k, m) / sumS;
}
```

**What it does:** Divides each bucket by that spectrum's total bucket sum, then scales by a hardcoded constant `100000.0`. Each spectrum is independently normalized to a fixed total of 100,000.

**Mathematical form:** `bucket[k,m] = 100000 × bucket[k,m] / sum(bucket[k,])`

### Why they differ

| Property | Path A (R) | Path B (C++) |
|---|---|---|
| Reference frame | Group mean | Per-spectrum sum |
| Relative magnitudes preserved across spectra? | Yes | No — every spectrum sums to 100,000 |
| Arbitrary scaling constant? | No | Yes (`100000.0`) |
| Operates on | Full spectral points | Integrated bucket values |

A user who expects "CSN" to mean the same thing in both places will get inconsistent results. The Path B approach (fixed-sum) is the more standard Constant Sum Normalization in metabolomics, but Path A's relative-to-mean approach preserves inter-spectrum magnitude differences, which is a meaningfully different operation.

### Recommendation

- Decide which CSN definition the project should use and implement it identically in both paths.
- If both approaches are intentionally kept, rename them (e.g., "CSN-relative" vs. "CSN-total") and document the difference for users.
- The `100000.0` magic constant in the C++ code should be configurable or at minimum documented.

---

## NORM-03: PQN Missing Required CSN Pre-Normalization (Path A)

| | |
|---|---|
| **Severity** | Critical |
| **Files** | `RnmrTools.R:437–451` |
| **Impact** | Spectral-level PQN is mathematically incorrect per the standard algorithm |
| **Reference** | Dieterle et al. (2006), Anal. Chem. 78(13):4281–4290 |

### The Standard PQN Algorithm

The Probabilistic Quotient Normalization algorithm as published by Dieterle et al. has four steps:

1. **Perform an initial integral normalization** (typically CSN) on all spectra
2. Calculate a reference spectrum (median of the CSN-normalized spectra)
3. Compute the quotient of each spectrum against the reference
4. Divide each spectrum by the median of its quotient vector

### What Path A does

```r
# RnmrTools.R, lines 437-451
if (normmeth=='PQN') { # TOTO: cf https://github.com/tkimhofer/metabom8/blob/master/R/pqn.R
    SUBMAT <- foreach(i=1:N, .combine=cbind) %dopar% {
        i1 <- length(which(specMat$ppm > max(zones[i,])))
        i2 <- which(specMat$ppm <= min(zones[i,]))[1]
        t(simplify2array(lapply(1:specMat$nspec, function(x) {
            specMat$int[x, i1:i2]
        })))
    }
    V <- apply(SUBMAT, 2, median)       # reference spectrum
    SUBMAT <- SUBMAT[, V != 0]
    V <- V[V != 0]
    MQ <- t(t(SUBMAT) / V)             # quotients
    COEFF <- apply(MQ, 1, median)       # median quotient per spectrum
}
```

**Step 1 is missing.** The raw, un-normalized spectra are used directly to compute the reference spectrum and quotients. Without the initial CSN step, spectra with very different total concentrations dominate the quotient calculation, making the PQN dilution factor unreliable.

Note the `# TOTO:` comment (likely meant `# TODO:`) referencing an external PQN implementation — this suggests the author was aware the implementation was incomplete.

### What Path B does (correctly)

```r
# RnmrTools.R, lines 1355-1360
if (norm_meth == 'PQN') {
    buckets_IntVal_CSN <- C_buckets_CSN_normalize(buckets_IntVal)  # Step 1 ✓
    bucVref_IntVal <- C_MedianSpec(buckets_IntVal_CSN)             # Step 2 ✓
    bucRatio <- sweep(buckets_IntVal_CSN, 2, bucVref_IntVal, "/")  # Step 3 ✓
    Coeff <- apply(bucRatio, 1, median)                            # Step 4 ✓
    buckets_IntVal <- sweep(buckets_IntVal_CSN, 1, Coeff, "/")
}
```

Path B correctly applies CSN first (`C_buckets_CSN_normalize`), then proceeds with steps 2–4. However, Path B has its own issue — it applies the final quotient division to the CSN-normalized data, not the original data. Standard PQN applies the median quotient to the original spectrum. This is a secondary concern but worth verifying against the reference implementation.

### Recommendation

Add the CSN pre-normalization step to Path A's PQN:

```r
if (normmeth=='PQN') {
    # Step 1: Initial CSN normalization
    # (compute SUM and COEFF as in the CSN branch above)
    SUM <- foreach(i=1:N, .combine='+') %dopar% { ... }
    CSN_COEFF <- SUM / mean(SUM)
    CSN_SUBMAT <- specMat$int / CSN_COEFF

    # Step 2-4: Standard PQN on the CSN-normalized spectra
    SUBMAT <- ...extract zones from CSN_SUBMAT...
    V <- apply(SUBMAT, 2, median)
    ...
    COEFF <- CSN_COEFF * apply(MQ, 1, median)
}
```

---

## NORM-04: Biased Median in `C_MedianSpec`

| | |
|---|---|
| **Severity** | High |
| **File** | `libCspec.cpp:457–470` |
| **Impact** | Systematic upward bias in PQN reference spectrum for even sample counts |

### Problem

```cpp
// libCspec.cpp, lines 462-468
int position = n_specs / 2;  // integer division
NumericVector out(count_max);
for (int j = 0; j < count_max; j++) {
    NumericVector y = VV(_, j);
    std::nth_element(y.begin(), y.begin() + position, y.end());
    out[j] = y[position];
}
```

For an even number of spectra `n`, the true median is `(element[n/2 - 1] + element[n/2]) / 2`. This code always takes `element[n/2]` — the upper of the two middle elements.

**Example:** For `n_specs = 10` and sorted values `{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}`:
- True median: `(5 + 6) / 2 = 5.5`
- This code returns: `6` (position = 5, 0-indexed)

This produces a systematically higher reference spectrum, which in turn makes all PQN quotients systematically smaller, which makes all PQN correction coefficients smaller, which makes all final normalized values **systematically larger** than they should be.

For typical NMR metabolomics datasets with 10–30 samples, this bias is meaningful and non-random.

### Recommendation

```cpp
int position = n_specs / 2;
if (n_specs % 2 == 0) {
    std::nth_element(y.begin(), y.begin() + position - 1, y.end());
    double lower = y[position - 1];
    std::nth_element(y.begin(), y.begin() + position, y.end());
    double upper = y[position];
    out[j] = (lower + upper) / 2.0;
} else {
    std::nth_element(y.begin(), y.begin() + position, y.end());
    out[j] = y[position];
}
```

---

## NORM-05: Division-by-Zero — No Guards Anywhere

| | |
|---|---|
| **Severity** | Critical |
| **Files** | Multiple (see below) |
| **Impact** | `NaN` / `Inf` values propagate silently into exported data matrices |

Every normalization path divides by a computed value. None check for zero.

### Location 1: `RNorm1D`, CSN — `RnmrTools.R:435`

```r
COEFF <- SUM / mean(SUM)
```

If all selected zones contain only zeros (e.g., user selects an empty region), `mean(SUM) = 0`, producing `NaN`. Then line 454 `MatInt <- specMat$int / COEFF` propagates `NaN` across the entire spectral matrix.

### Location 2: `RNorm1D`, PQN — `RnmrTools.R:449–450`

```r
MQ <- t(t(SUBMAT) / V)
COEFF <- apply(MQ, 1, median)
```

While `V != 0` is filtered on line 447–448, if ALL column medians are zero, `V` becomes empty, `SUBMAT` becomes a 0-column matrix, and `apply(MQ, 1, median)` returns `NaN` for every spectrum. There is no check for this edge case.

Additionally, if `COEFF` contains any zero values, line 454 `specMat$int / COEFF` produces `Inf`.

### Location 3: `C_buckets_CSN_normalize` — `libCspec.cpp:1088–1089`

```cpp
for (m = 0; m < n_bucs; m++) sumS += buckets(k, m);
for (m = 0; m < n_bucs; m++) M(k, m) = 100000.0 * buckets(k, m) / sumS;
```

If a spectrum's buckets sum to zero, `sumS = 0.0` and the division produces `Inf` or `NaN` (IEEE 754 `0.0/0.0 = NaN`). No check exists.

### Location 4: `get_Buckets_dataset`, PQN — `RnmrTools.R:1358`

```r
bucRatio <- sweep(buckets_IntVal_CSN, 2, bucVref_IntVal, "/")
```

If any element of the median reference spectrum `bucVref_IntVal` is zero, the sweep produces `Inf` in that column for all spectra. The `C_MedianSpec` function has no zero-filtering.

### Location 5: `get_Buckets_dataset`, Reference signal — `RnmrTools.R:1369`

```r
buckets_IntVal <- buckets_IntVal / Vref
```

If the reference signal integration `Vref` is zero for any spectrum (e.g., the selected reference region is outside the spectral range or contains only noise), that spectrum's bucket values become `Inf`.

### Recommendation

Add guards at every division point. Example pattern:

```r
if (any(COEFF == 0) || any(is.nan(COEFF)) || any(is.infinite(COEFF))) {
    bad_idx <- which(COEFF == 0 | is.nan(COEFF) | is.infinite(COEFF))
    Write.LOG(LOGFILE, paste0("WARNING: Normalization coefficient is zero/NaN/Inf for spectra: ",
                               paste(bad_idx, collapse=", ")))
    COEFF[bad_idx] <- 1.0  # leave these spectra unnormalized
}
```

For the C++ path, add a check before the division:

```cpp
if (sumS == 0.0) {
    for (m = 0; m < n_bucs; m++) M(k, m) = 0.0;  // or NaN, or flag
} else {
    for (m = 0; m < n_bucs; m++) M(k, m) = 100000.0 * buckets(k, m) / sumS;
}
```

---

## NORM-06: Silent Triple-Layer Normalization Stacking

| | |
|---|---|
| **Severity** | High |
| **Files** | `RnmrTools.R:1352–1369`, `ui_procparams.R:351–354` |
| **Impact** | Up to three normalization layers can be stacked without user awareness |

### The three layers

In `get_Buckets_dataset` (the export path), the following can all be active simultaneously:

**Layer 1** — The spectra may already be normalized from preprocessing (Path A via `RNorm1D`). The export function reads the already-modified `specs.pack`, so any prior spectral normalization is baked into the data.

**Layer 2** — Bucket-level CSN or PQN is applied on lines 1352–1361:
```r
if (norm_meth == 'CSN') {
    buckets_IntVal <- C_buckets_CSN_normalize(buckets_IntVal)
}
if (norm_meth == 'PQN') {
    ...
}
```

**Layer 3** — Reference signal normalization is applied on lines 1365–1369:
```r
if (! is.na(zoneref)) {
    Vref <- C_spectra_integrate(specMat$int, istart, iend)
    buckets_IntVal <- buckets_IntVal / Vref
}
```

All three layers are independent — there is no logic that checks whether one has been applied before applying another. The UI presents them as separate, unrelated options (the normalization dropdown and the "PPM range of the Reference" text box are in the same export panel).

### Example worst case

1. User applies PQN during preprocessing → spectra are PQN-normalized
2. User selects CSN in the export dropdown → bucket values are CSN-normalized on top of PQN
3. User enters a reference PPM range → a third normalization divides by the reference integral

The resulting values have been divided three times by three different factors, with no record of this in the exported data.

### Recommendation

- Expose the normalization status clearly in the export panel (e.g., "Spectral normalization applied: PQN").
- Disable or warn when stacking is detected.
- Consider a single normalization control point with clear "applied"/"not applied" state.

---

## NORM-07: Parallelization Bug in `RNorm1D`

| | |
|---|---|
| **Severity** | Medium |
| **File** | `RnmrTools.R:423–428` |
| **Impact** | Potential crash or undefined behavior with single-zone normalization |

### Problem

```r
RNorm1D <- function(specMat, normmeth, zones) {
    N <- dim(zones)[1]
    if (N > 1) registerDoParallel(cores = 2)  # Only registers when N > 1

    if (normmeth == 'CSN') {
        SUM <- foreach(i = 1:N, .combine = '+') %dopar% {  # Always uses %dopar%
            ...
        }
    }
}
```

When `N == 1` (single zone — a common case), `registerDoParallel` is not called. But the `foreach %dopar%` loop still expects a parallel backend. The behavior depends on whether a backend was previously registered elsewhere in the session:
- If yes: it works, but uses whatever core count was set before (potentially more than intended)
- If no: `%dopar%` may fall back to sequential or throw a warning depending on the `doParallel` version
- In either case: the behavior is **non-deterministic** across sessions

### Recommendation

Either always register the backend, or use `%do%` (sequential) when `N == 1`:

```r
if (N > 1) {
    registerDoParallel(cores = 2)
    SUM <- foreach(i = 1:N, .combine = '+') %dopar% { ... }
} else {
    SUM <- foreach(i = 1:N, .combine = '+') %do% { ... }
}
```

---

## NORM-08: Hardcoded Macro Replay Forces CSN

| | |
|---|---|
| **Severity** | Medium |
| **File** | `RnmrTools.R:1050–1059` |
| **Impact** | Macro replay ignores user-selected normalization method in legacy format |

### Problem

```r
if (cmdName == lbNORM) {
    params <- cmdPars[-1]
    if (length(params) == 2) {                                # Legacy 2-parameter format
        params <- as.numeric(params)
        PPMRANGE <- c(min(params[1:2]), max(params[1:2]))
        specMat <- RNorm1D(specMat, normmeth = 'CSN', ...)   # ← hardcoded 'CSN'
    }
    if (length(params) == 1) {                                # New 1-parameter format
        NORM_METH <- params[1]                                # ← reads method from macro
        specMat <- RNorm1D(specMat, normmeth = NORM_METH, ...)
    }
}
```

When replaying a macro command with the legacy 2-parameter format (`normalisation <ppm_min> <ppm_max>`), the method is hardcoded to `'CSN'` regardless of what the user originally selected. Only the newer 1-parameter format (`normalisation <method>` followed by zone lines) respects the saved method.

This means old macro files always replay as CSN, even if PQN was originally used.

### Recommendation

Either deprecate the legacy format with a clear warning, or add the method to the legacy format.

---

## NORM-09: Inconsistent Parameter Naming

| | |
|---|---|
| **Severity** | Low (maintainability) |
| **Files** | Multiple |
| **Impact** | Developer confusion; risk of wiring errors |

The normalization method parameter appears under three different names:

| Context | Variable Name | File |
|---|---|---|
| Preprocessing UI | `input$normeth` | `ui_procparams.R:129` |
| Export UI | `input$normmeth` | `ui_procparams.R:351` |
| Processing params | `procParams$NORM_METH` | `Proc4.R:23` |
| Function argument | `normmeth` | `RnmrTools.R:421` |
| Function argument | `norm_meth` | `RnmrTools.R:1329` |

This makes it difficult to trace which normalization setting is being used where, and invites copy-paste bugs.

### Recommendation

Standardize on a single name (e.g., `norm_method`) across the entire codebase.

---

## General Issues

### GEN-01: `eval(parse(...))` Throughout INI Parsing — Security Risk

**File:** `RnmrTools.R:48, 79, 83, 86`

The `Parse.INI` and `Write.INI` functions construct R expressions from file contents and execute them:

```r
eval(parse(text = paste0('INI.list$', d$V1[i], '<-"', d$V2[i], '"')))
```

Since the INI file is generated from user input (via `generate_INI_file` in `Proc4.R`), a crafted value could inject arbitrary R code. Safe replacement:

```r
INI.list[[d$V1[i]]] <- d$V2[i]
```

### GEN-02: Shell Injection in 7z Extraction

**File:** `RnmrTools.R:153`

```r
system(paste0("cd ", RAWDIR, "; 7zr x -y ", RawZip))
```

`RawZip` derives from the uploaded filename. A file named `; rm -rf /; .7z` would execute arbitrary commands. Use `system2("7zr", args = c("x", "-y", RawZip), ...)` instead.

### GEN-03: No Unit Tests

There are zero test files in the repository. For numerical processing code where bugs silently corrupt results (as demonstrated by the normalization issues), automated testing is essential. At minimum, normalization functions should have tests for:
- Known-answer cases (e.g., all-equal spectra should be unchanged)
- Edge cases (single spectrum, single bucket, zero values, negative values)
- Consistency between Path A and Path B for the same input
- Regression cases from reported bugs

### GEN-04: Monolithic `RnmrTools.R` (1,459 Lines)

This single file contains: INI parsing, logging, stack management, metadata generation, baseline correction (3 methods), calibration, normalization, alignment, bucketing, filtering, smoothing, spectra I/O, and bucket dataset generation. It should be broken into focused modules for maintainability.

### GEN-05: Duplicated Bucket Parsing Logic

The bucket file (`bucket_list.in`) is read and the `min/max` from `center/width` is recomputed in at least four places:
- `get_Buckets_table` (line 1309)
- `get_Buckets_dataset` (line 1334)
- `get_SNR_dataset` (line 1418)
- `Rbuc1D` (lines 412, 465)

This should be a single shared function.

### GEN-06: Warnings Suppressed Globally

**File:** `Rcorr1D:14` — `options(warn = -1)`

This suppresses all R warnings for the entire preprocessing session. Any numerical issues (NaN produced, coercion warnings, convergence warnings) are silently swallowed. This directly masks the division-by-zero issues in NORM-05.

### GEN-07: Hardcoded Noise Range Default

**Files:** `Rcorr1D` (multiple locations), `Rbuc1D:129`

```r
PPM_NOISE_AREA <- c(10.2, 10.5)  # ¹H default
```

This default is appropriate for proton NMR but incorrect for ¹³C (where the config already defines `PPM_MIN_13C=-20, PPM_MAX_13C=200`). The noise range default should be nucleus-aware.

### GEN-08: Typo in UI

**File:** `ui_procparams.R:122–123`

```r
numericInput("ppmshift2", "PPM shitf value:", ...)
bsTooltip("ppmshift2", "PPM shitf value to be applied ...", ...)
```

"shitf" → "shift" (appears twice).

---

## Proposed Fix Priority

### Tier 1 — Data Correctness (fix before any release)

| ID | Issue | Effort |
|---|---|---|
| NORM-03 | Add CSN pre-normalization to Path A PQN | Medium |
| NORM-05 | Add zero-division guards to all 5 normalization division points | Medium |
| NORM-04 | Fix `C_MedianSpec` even-count median calculation | Low |
| NORM-02 | Unify CSN algorithm between Path A and Path B (or rename/document) | Medium |
| GEN-06 | Remove `options(warn=-1)` — at minimum, log warnings to file | Low |

### Tier 2 — User Protection (fix in next iteration)

| ID | Issue | Effort |
|---|---|---|
| NORM-01 | Prevent or warn about double normalization | Medium |
| NORM-06 | Add normalization status tracking across layers | Medium |
| NORM-08 | Fix hardcoded CSN in macro replay | Low |
| NORM-09 | Standardize parameter naming | Low |

### Tier 3 — Code Health (ongoing improvement)

| ID | Issue | Effort |
|---|---|---|
| GEN-03 | Add unit tests for normalization functions | High |
| GEN-01 | Replace `eval(parse(...))` with safe `[[` assignment | Low |
| GEN-02 | Sanitize shell commands | Low |
| GEN-04 | Refactor `RnmrTools.R` into modules | High |
| GEN-05 | Extract shared bucket parsing | Low |
| NORM-07 | Fix parallelization for single-zone case | Low |
| GEN-07 | Nucleus-aware noise defaults | Low |
| GEN-08 | Fix "shitf" typo | Trivial |

---

*End of review. Questions or clarifications: reach out to the reviewing team.*
