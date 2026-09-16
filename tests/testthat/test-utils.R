test_that("data trimming helpers", {
  harmonized <- harmonize_features(list(A=data.frame(t1=1:3, t2=4:6), B=data.frame(t2=1:2, t3=3:4)), verbose=FALSE)
  expect_identical(colnames(harmonized$B), "t2")
  md <- normalize_NA_in_metadata(data.frame(a=c("x", "not provided"), b=factor(c("NA", "y"))))
  expect_true(is.na(md$a[2]) && is.na(md$b[1]))
  expect_equal(ncol(discard_uninfo_columns_in_metadata(data.frame(a=1:2, b=3:4))), 2)
  expect_equal(dim(filter_features_allzero(data.frame(a=c(0, 0, 1), b=c(0, 0, 0), c=c(0, 1, 1)))), c(2, 2))
  expect_s3_class(check_metadata(data.frame(a=c("x", "?"), b=1:2), more_missing_values="?"), "data.frame")
})
