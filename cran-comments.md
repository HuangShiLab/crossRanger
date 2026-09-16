# cran-comments

## Submission

This is a new submission of crossRanger to CRAN.

crossRanger builds random forest models (using 'ranger') within multiple microbiome datasets,
transfers the models across datasets, and identifies biomarkers that generalize across studies.

## Test environments

* local: macOS 26.5.1 (aarch64-apple-darwin20), R 4.3.3

Before submitting, please add the results of at least:

* win-builder (R-release and R-devel)
* R-hub: Windows Server, Ubuntu Linux and Fedora

## R CMD check results

`R CMD check --as-cran` gives:

```
0 errors | 0 warnings | 2 notes
```

### Note 1: checking CRAN incoming feasibility

```
Maintainer: 'Shi Huang <shihuang047@gmail.com>'

New submission
```

This is expected for a first submission.

The same note also reported one URL as possibly invalid:

```
URL: https://www.gnu.org/licenses/gpl-3.0.en.html
  From: README.md
  Status: Error
  Message: libcurl error code 35:
    LibreSSL SSL_connect: SSL_ERROR_SYSCALL in connection to www.gnu.org:443
```

This is a TLS failure of the machine the check was run on rather than a broken link: `curl`
fails with the same error for this address on that machine, while other https addresses (for
example <https://cran.r-project.org>) return HTTP 200 from the same session. The address is the
canonical page of the GPL-3 license and resolves normally elsewhere.

### Note 2: checking for future file timestamps

```
unable to verify current time
```

The machine running the check could not reach the time server used by `R CMD check`. This is a
property of the check environment and unrelated to the package.

## Notes on the examples

* The examples of `crossRanger_report()` are wrapped in `\dontrun{}`: rendering an R Markdown
  report requires pandoc and writes an output file.
* A few examples that fit many random forest models are wrapped in `\donttest{}` to keep the
  check time short. They are run and pass with `R CMD check --run-donttest`.
* No example, test or vignette writes outside `tempdir()`.

## Downstream dependencies

There are currently no downstream dependencies (new submission).
