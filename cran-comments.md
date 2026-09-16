# cran-comments

## Submission

This is a new submission of crossRanger to CRAN.

crossRanger builds random forest models (using 'ranger') within multiple microbiome datasets,
transfers the models across datasets, and identifies biomarkers that generalize across studies.

## Test environments

* win-builder: R-devel (2026-09-15 r90540 ucrt), x86_64-w64-mingw32
* local: macOS 26.5.1 (aarch64-apple-darwin20), R 4.3.3

## R CMD check results

On win-builder (R-devel) the check gives one note:

```
Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Shi Huang <shihuang047@gmail.com>'

  New submission
```

This note is expected for a first submission; there are no errors or warnings.

Locally (macOS, R 4.3.3) `R CMD check` gives `Status: OK`.

## Notes on the examples

* The examples of `crossRanger_report()` are wrapped in `\dontrun{}`: rendering an R Markdown
  report requires pandoc and writes an output file.
* A few examples that fit many random forest models are wrapped in `\donttest{}` to keep the
  check time short. They are run and pass with `R CMD check --run-donttest`.
* No example, test or vignette writes outside `tempdir()`.

## Parallel computation

The package parallelizes over datasets with foreach and doParallel, which forks on Unix and uses
a PSOCK cluster on Windows. A cluster started by the package is stopped again by the function that
started it, so that no connection is left open, and a parallel backend registered by the user is
used as it is and never stopped. While the package is checked the functions run sequentially.

## Downstream dependencies

There are currently no downstream dependencies (new submission).
