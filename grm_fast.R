#' Calculates the genomic relationship matricx using the first method of VanRaden 2008.
#' Expects a dataframe or datatable with animal identificaton as the first column.
#' CANNOT yet handle missing values or filtering.
#'
#' @param M dataframe or datatable with marker information and animal identification
#' @return Genomic relationship matrix
#' @export
#'

grm.fast = function(M, ids = T, groups = F){
  if (!require('data.table', quietly = TRUE)) { install.packages('data.table') }
  library('data.table')

  if(!("data.table") %in% class(M)){
    M = as.data.table(M)
  }

  if(ids == T && groups == T){
    groups_and_id = as.data.frame(M[,1:2])
    M[,colnames(M)[1:2]:=NULL]
  }
  if(ids == T && groups == F){
    saved_names = as.data.frame(M[,1])
    M[,colnames(M)[1]:=NULL]
  }

  if(("matrix") %in% class(M)){
    cat("Object is matrix. \n")
  }else{
    M = as.matrix(M)
    cat("Object converted to matrix. \n")
  }

  p = colSums(M)/(2*nrow(M))
  k = 2*sum(p*(1-p))

  P = matrix(rep(p, each = nrow(M)),
             ncol = ncol(M))
  P = 2*P
  Z = M - P
  rm(M, P, p)

  G = tcrossprod(Z)/k
  colnames(G) = rownames(G) = saved_names[,1]

  return(G)
}
