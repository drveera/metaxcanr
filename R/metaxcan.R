#' metaxcan main function
#'
#' @param gwas.file gwas summary stats file name
#' @param db.file prediction model db file name
#' @param snpcov.file pairwise snp-covariance file name
#'
#' @import DBI
#' @import dbplyr
#' @import dplyr
#' @import RSQLite
#' @import data.table
#' @import foreach
#' @import doParallel
#' @import parallel
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

  ##load the db file
  db.con <- dbConnect(SQLite(), db.file)
  db.df <- tbl(db.con,"weights")  %>% collect()
  db.df <- data.table(db.df)
  genes <- unique(db.df$gene)
  res <- list()
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  res <- foreach(i=genes,
                 .combine = rbind,
                 .packages=c("data.table","metaxcanr")) %dopar%
    impute.zscore(gene.name=i,gwas=gwas,db=db.df,snpcov=snpcov)
  ##get extras
  extra <- tbl(db.con,"extra") %>% collect()
  extra <- data.table(extra)
  res <- merge(res,extra,by="genename")
  return(res)
}

##library(DBI,lib.loc = "~/va-biobank/Veera/Rlibraries/")
##library(dbplyr)
##library(dbplyr,,lib.loc = "~/va-biobank/Veera/Rlibraries/")
##library(dplyr,lib.loc = "~/va-biobank/Veera/Rlibraries/")
##library(RSQLite, lib.loc = "~/va-biobank/Veera/Rlibraries/")
##library(data.table)
##library(foreach)
##library(doParallel)
