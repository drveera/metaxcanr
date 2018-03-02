#' impute the z score
#'
#'
#' @param gene.name gene name
#' @param gwas gwas data frame
#' @param db data frame containing the prediction models
#' @param snpcov data frame containing the pairwise snp covariances
#'
#'
#' @export
impute.zscore <- function(gene.name,gwas,db,snpcov){
  ##subset gene information
  ##db <- tbl(db,'weights') %>% filter(gene==gene.name) %>% collect()
  db <- db[gene == gene.name]
  n_snps_in_model <- nrow(db)
  snpcov <- snpcov[GENE==gene.name]
  gwas <- gwas[SNP %in% db$rsid]
  n_snps_used <- nrow(gwas)
  if(nrow(gwas)<=1){
    zscore=NA
    effsize=NA
    res <- data.frame(gene=gene.name,
                      zscore=zscore,
                      effsize=effsize,
                      n_snps_used=n_snps_used,
                      n_snps_in_model=n_snps_in_model,
                      stringsAsFactors = FALSE)
    return(res)
  }

  ##cov matrix
  snpcov.mat <- create.snpcov.matrix(snpcov)
  snps <- rownames(snpcov.mat)
  snps <- snps[snps %in% gwas$SNP]
  ##update matrix
  snpcov.mat <- snpcov.mat[snps,snps]
  ##keep weights only present in gwas and cov mat
  db <- db[rsid %in% snps]

  ##cat("total no. of snps used for",gene.name,nrow(gwas),"\n")

  ##weights
  snpwts <- db$weight
  names(snpwts) <- db$rsid
  snpwts <- snpwts[snps]

  ##effalleles
  effalleles <- db$eff_allele
  names(effalleles) <- db$rsid
  effalleles <- effalleles[snps]

  ##betas
  betas <- gwas$BETA
  names(betas) <- gwas$SNP
  betas <- betas[snps]

  ##se
  serrors <- gwas$SE
  names(serrors) <- gwas$SNP
  serrors <- serrors[snps]

  ##zscores
  zscores <- betas/serrors
  names(zscores) <- snps

  ##a1
  a1 <- gwas$A1
  names(a1) <- gwas$SNP
  a1 <- a1[snps]

  ##update zscores and betas
  zscores <- ifelse(a1==effalleles,zscores,(zscores*-1))
  betas <- zscores*serrors


  ##get diagonals
  snpcov.diag <- diag(snpcov.mat)


  ##sigmas
  sigmas <- sqrt(snpcov.diag)

  ##calculate
  cat(length(snpwts),dim(snpcov.mat),"\n")
  genevariance <- (snpwts %*% snpcov.mat) %*% snpwts

  zscore <- sum(snpwts * zscores * sigmas)/sqrt(genevariance)
  effsize <- sum(snpwts * betas * (sigmas^2))/genevariance

  res <- data.frame(gene=gene.name,
                    zscore=zscore,
                    effsize=effsize,
                    n_snps_used=n_snps_used,
                    n_snps_in_model=n_snps_in_model,
                    stringsAsFactors = FALSE)
  return(res)
}
