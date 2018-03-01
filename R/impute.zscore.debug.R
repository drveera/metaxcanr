#' debug mode for impute zcore function
#'
#'
#'
#'
#' @export
impute.zscore.debug <- function(geneid,...){
  tryCatch({
    impute.zscore(...)
  },
  error = function(e){
    return(data.frame(gene=geneid,
           zscore=NA,
           effsize=NA,
           n_snps_used=NA,
           n_snps_in_model=NA,
           error_message = paste(e),
           stringsAsFactors = FALSE))
  })
}
