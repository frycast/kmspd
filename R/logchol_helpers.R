library(abind)
library(reticulate)

## Map data to the Cholesky space
# Convention here is that the first dim*(dim-1)/2 columns are the triangle
# and the remaining dim columns are the diagonal. The attribute tr will
# indicate the columns that belong to lower triangles.
map_to_cholesky_log <- function(data) {
  lower_tr_mat <- lower.tri(data[,,1])
  d <- dim(data[,,1])[1]
  n <- dim(data[1,,])[2]
  dataL <- array(dim = dim(data))
  for (i in 1:n) {dataL[,,i] <- t(chol(data[,,i]))}
  dataL_tr <- t(apply(dataL, FUN = function(x){x[lower_tr_mat]}, MARGIN = 3))
  dataL_di <- log(t(apply(dataL, FUN = diag, MARGIN = 3)))
  dataL_flat <- cbind(dataL_tr, dataL_di)
  #attr(dataL_flat, "lower_tr") <- c(rep(T, sum(lower_tr_mat)), rep(F, d))
  return(dataL_flat)
}

# New version of the above that works with list of matrices instead of array
map_to_cholesky_log2 <- function(data) {
  lower_tr_mat <- lower.tri(data[[1]])
  d <- nrow(data[[1]])
  n <- length(data)
  dataL <- list()
  dataL <- lapply(data, function(x){t(chol(x))}) 
  dataL_flat <- t(sapply(dataL, 
    FUN = function(x){
      tr <- x[lower_tr_mat]
      di <- log(diag(x))
      c(tr, di)
    }, simplify = "array"))
  return(dataL_flat)
}

map_from_cholesky_log <- function(data) {
  n <- nrow(data)
  nc <- ncol(data)
  d <- (sqrt(1+8*nc) - 1)/2
  lower_tr <-  c(rep(T, choose(d,2)), rep(F, d))
  data[,!lower_tr] <- exp(data[,!lower_tr])
  data_mat <- array(0L, dim = c(d,d,n))
  lt <- lower.tri(data_mat[,,1])
  for (k in 1:n) {
    diag(data_mat[,,k]) <- data[k,!lower_tr]
    data_mat[,,k][lt] <- data[k,lower_tr]
    data_mat[,,k] <- data_mat[,,k] %*% t(data_mat[,,k])
  }
  return(data_mat)
}