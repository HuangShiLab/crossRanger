#' @title crossRanger_report
#' @description Generate a dynamic, reproducible R Markdown report containing the analysis results,
#' performance summaries and figures of crossRanger analyses, rendered to HTML, PDF or Word.
#' @param ... One or more (optionally named) result objects, e.g., outputs of \code{rf_clf.by_datasets},
#' \code{rf_reg.by_datasets}, \code{rf_clf.by_datasets.summ}, \code{rf_clf.comps.summ}, \code{rf_clf.cross_appl},
#' \code{rf_reg.cross_appl}, \code{rf_clf.lodo}, \code{rf_reg.lodo}, \code{rf_clf.rfe}, \code{rf_reg.rfe},
#' \code{corr_datasets_by_imps}, a \code{BetweenGroup.test} data.frame, or ggplot objects.
#' The names are used as section titles.
#' @param output_format The R Markdown output format: "html_document" (default), "pdf_document" (requires LaTeX)
#' or "word_document".
#' @param output_file The file name of the report. By default "crossRanger_report" with the extension of the format.
#' @param outdir The output directory. By default, the working directory.
#' @param title The title of the report.
#' @param top_n The number of top-ranked features listed in the tables.
#' @return The path of the rendered report (invisibly).
#' @seealso rf_clf.by_datasets.summ rf_clf.comps.summ
#' @examples
#' \dontrun{
#' data(crossRanger_sim)
#' res <- rf_clf.by_datasets(crossRanger_sim$abundance, crossRanger_sim$metadata,
#'                           s_category="cohort", c_category="disease_status", positive_class="Case")
#' cross_rf <- rf_clf.cross_appl(res$rf_model_list, res$x_list, res$y_list, positive_class="Case")
#' crossRanger_report(Within_study_models=res, Cross_application=cross_rf, outdir=tempdir())
#' }
#' @export
crossRanger_report <- function(..., output_format="html_document", output_file=NULL, outdir=NULL,
                               title="crossRanger analysis report", top_n=20){
  if(!requireNamespace("rmarkdown", quietly=TRUE)) stop("Package 'rmarkdown' is required to render reports.")
  output_format <- match.arg(output_format, c("html_document", "pdf_document", "word_document"))
  results <- list(...)
  if(length(results)==0) stop("Provide at least one crossRanger result object.")
  result_names <- names(results)
  if(is.null(result_names)) result_names <- rep("", length(results))
  unnamed <- result_names==""
  result_names[unnamed] <- vapply(results[unnamed], function(r) class(r)[1], character(1))
  names(results) <- make.unique(result_names)
  extension <- c(html_document="html", pdf_document="pdf", word_document="docx")[[output_format]]
  if(is.null(output_file)) output_file <- paste0("crossRanger_report.", extension)
  if(is.null(outdir)) outdir <- getwd()
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  template <- system.file("rmarkdown", "crossRanger_report.Rmd", package="crossRanger")
  if(!nzchar(template)) stop("The report template was not found. Please reinstall crossRanger.")
  work_dir <- tempfile("crossRanger_report_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive=TRUE), add=TRUE)
  file.copy(template, work_dir)
  out <- rmarkdown::render(file.path(work_dir, basename(template)), output_format=output_format,
                           output_file=output_file, output_dir=normalizePath(outdir),
                           intermediates_dir=work_dir,
                           params=list(results=results, title=title, top_n=top_n),
                           envir=new.env(parent=globalenv()), quiet=TRUE)
  message("The report was saved to: ", out)
  invisible(out)
}
