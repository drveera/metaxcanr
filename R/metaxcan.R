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
                     genes=NA,
                     ncores=1){
  ##load the db file
  db.con <- dbConnect(SQLite(), db.file)
  db.df <- tbl(db.con,"weights")  %>% collect()
  db.df <- data.table(db.df)
  if(is.na(genes)[1]){
    genes <- unique(db.df$gene)
  } else {
    norig <- length(genes)
    genes <- intersect(unique(db.df$gene),genes)
    if(length(genes)==0){
      warning("None of the input genes present in db file. Nothing returned")
      return(NULL)
    }
    cat("Imputing ",length(genes),"genes out of ",norig," user provided genes \n")
    ##read the gwas
    gwas <- fread(gwas.file, header=TRUE)
    ##keep only SNP, A1, A2, BETA, SE
    gwas <- gwas[,c("SNP","A1","A2","BETA","SE"),with=FALSE]

    ##read the covariance file
    snpcov <- fread(snpcov.file, header=TRUE)
    snpcov <- snpcov[!is.na(VALUE)]


###
    ##subset gwas
    gwas <- gwas[SNP %in% db.df$rsid]
    ##keep only genes with predictors in gwas file
    db.df <- db.df[rsid %in% gwas$SNP]
    snpcov <- snpcov[GENE %in% db.df$gene]
####
    res <- list()
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    res <- foreach(i=genes,
                   .combine = bind_rows,
                   .packages=c("data.table","metaxcanr","dplyr")) %dopar%
      impute.zscore.debug(geneid=i,gene.name=i,gwas=gwas,db=db.df,snpcov=snpcov)
    res$pvalue <- 2 * pnorm(-abs(res$zscore))
  }
  ##get extras
  extra <- tbl(db.con,"extra") %>% collect()
  extra <- data.table(extra)
  res <- merge(res,extra,by="gene")
  dbDisconnect(db.con)
  return(res)
}
