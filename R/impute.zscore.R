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
impute.zscore <- function(gene.name,
                          gwas,db,
                          snpcov,
                          snpinfo,
                          zscore=NA){
  ##subset gene information
  ##db <- tbl(db,'weights') %>% filter(gene==gene.name) %>% collect()
  db <- db[gene == gene.name]
  n_snps_in_model <- nrow(db) ##delete
  snpcov <- snpcov[GENE==gene.name]
  gwas <- gwas[SNP %in% db$rsid]
  ##extract the minimum p value
  minP <- min(gwas$P,na.rm = TRUE)
  n_snps_used <- nrow(gwas) ##keep this for using it in NA genes
  if(nrow(gwas)<=1){
    zscore=NA
    effsize=NA
    res <- data.frame(gene=gene.name,
                      zscore=zscore,
                      effsize=effsize,
                      n_snps_used=n_snps_used,
                      n_snps_noCov_info = NA,
                      gwas_minP = minP,
                      stringsAsFactors = FALSE,
                      error=paste0("No. of GWAS SNPs=",nrow(gwas)))
    return(res)
  }
  ##cov matrix
  snpcov.mat <- create.snpcov.matrix(snpcov)
  snps <- rownames(snpcov.mat)
  snps <- snps[snps %in% gwas$SNP]
  n_snps_used <- length(snps) ##this is actual n
  n_snps_noCov_info <- n_snps_in_model - n_snps_used
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

  if(is.na(zscore)){
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
  } else {
    zscores <- gwas$Z
    names(zscores) <- gwas$SNP
    zscores <- zscores[snps]
  }

  ##a1
  a1 <- gwas$A1
  names(a1) <- gwas$SNP
  a1 <- a1[snps]

  ##update zscores and betas
  zscores <- ifelse(a1==effalleles,zscores,(zscores*-1))
  if(is.na(zscore)) betas <- zscores*serrors


  ##get diagonals
  snpcov.diag <- diag(snpcov.mat)


  ##sigmas
  sigmas <- sqrt(snpcov.diag)

  ##calculate
  cat(length(snpwts),dim(snpcov.mat),"\n")
  genevariance <- (snpwts %*% snpcov.mat) %*% snpwts

  zscore.gene <- sum(snpwts * zscores * sigmas)/sqrt(genevariance)
  if(is.na(zscore)){
    effsize <- sum(snpwts * betas * (sigmas^2))/genevariance
  } else {
    effsize=NA
  }

  res <- data.frame(gene=gene.name,
                    zscore=zscore.gene,
                    effsize=effsize,
                    n_snps_used=n_snps_used,
                    n_snps_noCov_info = n_snps_noCov_info,
                    gwas_minP = minP,
                    stringsAsFactors = FALSE)
  if(snpinfo) res$snps_used_for_prediction=paste(snps,collapse =";")
  return(res)
}
