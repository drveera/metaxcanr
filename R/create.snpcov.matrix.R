#' create a matrix from pairwise data
#'
#'
#' @param snpcov pairwise covariance matrix
#' @param snps a vector of snps to subset
#'
#' @export
create.snpcov.matrix <- function(snpcov,snps){
  snpcov.dup <- snpcov
  snpcov.dup$RSID1 <- snpcov$RSID2
  snpcov.dup$RSID2 <- snpcov$RSID1
  snpcov <- rbind(snpcov,snpcov.dup)
  ##to matrix
  snpcov.mat <- reshape(snpcov,idvar="RSID1",timevar = "RSID2", direction="wide")
  names(snpcov.mat) <- gsub("VALUE.","",names(snpcov.mat))
  snp.ids <- snpcov.mat$RSID1
  snpcov.mat <- snpcov.mat[,-1]
  snpcov.mat <- sapply(snpcov.mat,as.numeric)
  snpcov.mat <- as.matrix(snpcov.mat)
  rownames(snpcov.mat) <- snp.ids
  snpcov.mat <- snpcov.mat[snps,snps]
  return(snpcov.mat)
}
