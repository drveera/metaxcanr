#' plot weights of predictors and highlight a snp of interest
#'
#' @param snp snp id
#' @param gene gene id
#' @param db.file db file name with full path
#'
#' @importFrom DBI dbConnect dbDisconnect
#' @import dbplyr
#' @import dbplyr
#' @importFrom RSQLite SQLite
#' @import data.table
#'
#'@export
plot.wts <- function(snp,gene,db.file){
  ##load the db file
  db.con <- dbConnect(SQLite(), db.file)
  db.df <- tbl(db.con,"weights")  %>% collect()
  db.df <- data.table(db.df)
  geneid <- gene
  db.df <- db.df[gene==geneid]
  if(snp %in% db.df$rsid){
    message("Query SNP is not in the predictors")
    return(NULL)
  }
  snp.df <- db.df[rsid==snp]
  db.df <- db.df[order(weight,decreasing = TRUE)]
  db.df$rsid <- factor(db.df$rsid, levels=db.df$rsid)
  ##plot
  p1 <- ggplot(db.df, aes(rsid,weight)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_point(data=snp.df, fill="red",shape=23, size=2) +
    geom_text(data=snp.df,label=snp.df$rsid, angle=90, nudge_y=0.1)+
    ylab(snp.df$gene)
  return(p1)
}

