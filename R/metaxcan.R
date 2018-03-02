#' metaxcan main function
#'
#' @param gwas.file gwas summary stats file name
#' @param db.file prediction model db file name
#' @param snpcov.file pairwise snp-covariance file name
#'
#' @importFrom  DBI dbConnect dbDisconnect
#' @import dbplyr
#' @import dplyr
#' @importFrom RSQLite SQLite
#' @import data.table
#' @import foreach
#' @import doParallel
#' @importFrom parallel makeCluster
#' 
#' @export
metaxcan <- function(gwas.file,
                     db.file,
                     snpcov.file,
                     ncores=1){
  ##read the gwas
  gwas <- fread(gwas.file, header=TRUE)
  ##keep only SNP, A1, A2, BETA, SE
  gwas <- gwas[,c("SNP","A1","A2","BETA","SE"),with=FALSE]

  ##read the covariance file
  snpcov <- fread(snpcov.file, header=TRUE)
  snpcov <- snpcov[!is.na(VALUE)]

  ##load the db file
  db.con <- dbConnect(SQLite(), db.file)
  db.df <- tbl(db.con,"weights")  %>% collect()
  db.df <- data.table(db.df)

###
  ##subset gwas
  gwas <- gwas[SNP %in% db.df$rsid]
  ##keep only genes with predictors in gwas file
  db.df <- db.df[rsid %in% gwas$SNP]
  snpcov <- snpcov[GENE %in% db.df$gene]
####

  genes <- unique(db.df$gene)
  res <- list()
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  res <- foreach(i=genes,
                 .combine = bind_rows,
                 .packages=c("data.table","metaxcanr","dplyr")) %dopar%
    impute.zscore.debug(geneid=i,gene.name=i,gwas=gwas,db=db.df,snpcov=snpcov)
  ##get extras
  extra <- tbl(db.con,"extra") %>% collect()
  extra <- data.table(extra)
  res$pvalue <- 2 * pnorm(-abs(res$zscore))
  res <- merge(res,extra,by="gene")
  dbDisconnect(db.con)
  return(res)
}
