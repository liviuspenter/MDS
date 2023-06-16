# source:
# https://github.com/caleblareau/asap_reproducibility/blob/master/pbmc_stim_multiome/code/11_setup.R

import_kite_counts2 <- function(dir){
  mtx <- fread(paste0(dir, "/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0(dir, "/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0(dir, "/featurecounts.genes.txt"), header = FALSE)[[1]])[1:ncol(matx)]
  return(t(matx))
}
