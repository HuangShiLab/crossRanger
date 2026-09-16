.onUnload <- function(libpath){
  # On Windows, doParallel::registerDoParallel(cores=) starts an implicit PSOCK cluster,
  # whose worker processes have to be stopped explicitly.
  if(requireNamespace("doParallel", quietly=TRUE)) try(doParallel::stopImplicitCluster(), silent=TRUE)
}
