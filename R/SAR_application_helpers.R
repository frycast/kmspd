###########################################################################
# MAIN INTERFACE ----------------------------------------------------------
###########################################################################

#' get_kms_grid
#'
#' Run k-means for the grid (BAND, D, M, k), for each k in {1,2,...,k_max},
#' saving all results along the way.
#' 
#' On subsequent runs the loop will just load in results,
#' where the results are stored in the SAR_app_data folder. 
#' 
#' @return From the
#' results, using result_fun, are extracted, 
#' and, if format = "df", a data.frame is returned that is suitable
#' for use with ggplot2. Otherwise, if format = "list", a 
#' list of ssq vectors are returned,
#' where the length of each vector is k_max and the list element names
#' identify the triple (band, d, m).
#' Note the entire kms object is saved to disk, 
#' along with intermediate data at many steps (e.g., autocovariance matrices),
#' so that it can all be accessed later if needed.
#' 
#' @param BAND character vector of band names, from "vv", "vh" and "cc"
#' @param D integer vector of patch sizes
#' @param M integer vector of lags
#' @param k_max maximum number of clusters
#' @param iter.max argument to kmeans {stats}
#' @param nstart argument to kmeans {stats}
#' @param data_file the data file specifically for this function
#' @param region the region of interest (as appearing in the tif file name)
#' @param verbose whether to output status reports during execution
#' @param save whether to save results that are not already saved
#' @param result_fun a function that will be applied, using sapply, to kms.
#' result_fun should return one atomic type (e.g., numeric) for each
#' value of k in each grid point.
#' The default result_fun returns the total within cluster sum of squares. 
#' Labels can be incorporated as the second argument to result_fun if
#' labs_in_rf is set to TRUE.
#' @param format the output format, either "df" or "list"
#' @param labs_in_rf whether to load labels and use them in 
#' result_fun at each grid point 
get_kms_result <- function(BAND = c("cc","vv","vh"),
                         D = c(4,6,7,8,9,10),
                         M = c(1,2,3,4,5),
                         k_max = 50, iter.max = 100, nstart = 20,
                         data_file = "SAR_app_data", region = "MG",
                         verbose = T, save = T,
                         result_fun = function(kms){kms$tot.withinss},
                         format = "df",
                         labs_in_rf = F) {
  res <- list(); i <- 1
  for (d in D) {
    if (labs_in_rf) {
      labs <- get_pix_patched("labs", d, save, data_file, region, verbose)
    }
    for (band in BAND) {
      for (m in M) {
        if (verbose) {
          this_set <- paste0(band,": ","d = ",d,", m = ",m)
          message("-------------------------")
          message(this_set)
          message("-------------------------")
        }
        kms <- get_kms(band, d, m, save, data_file, region, 
                       k_max, nstart, iter.max, verbose)
        if (labs_in_rf) {
          res[[i]] <- sapply(kms, result_fun, labs)
        } else {
          res[[i]] <- sapply(kms, result_fun)
        }
        i <- i + 1
      }
    }
  }
  grid_labels <- expand.grid(M,BAND,D)
  lnams <- apply(grid_labels, function(x){
    paste0("m",x[1],"b",x[2],"d",x[3])}, MARGIN = 1)
  lnams <- gsub(" ", "0", lnams)
  if (format == "df") {
    
    grid_labels <- expand.grid(1:k_max,M,BAND,D)
    res <- do.call(c, res)
    df <- cbind(res, grid_labels, rep(lnams, each = k_max))
    colnames(df) <- c("res","k","m","band","d","run")
    return(df)
  } else if (format == "list") {
    
    names(res) <- lnams 
    return(res)
  } else {
    stop("format not recognised") 
  }
}

#' get_kms
#'
#' Given pixels, labels, patch size and lag, 
#' repeat k-means with penalty used to decide on 
#' optimal k up to k_max. 
#' 
#' @param band the band, vv, vh or cc
#' @param d patch size
#' @param m lag
#' @param save whether to save results that are not already saved
#' @param data_file the data file specifically for this function
#' @param region the region of interest (as appearing in the tif file name)
#' @param k_max maximum number of clusters
#' @param nstart argument to kmeans {stats}
#' @param iter.max argument to kmeans {stats}
#' @param verbose whether to output status reports during execution
get_kms <- function(band, d, m, save = T, data_file = "SAR_app_data", 
                    region = "MG", 
                    k_max = 50, nstart = 20, 
                    iter.max = 100, verbose = T) {
  
  fn <- paste0(data_file,"/",toupper(region),"/",
               "/kms_d",d,"_m",m,"_",band,
               "_kmax",k_max,"_ns",nstart,"_mi",iter.max,".Rds")
  if (file.exists(fn)) {
    set <- paste0("band = ",band,", d = ",d, ", m = ",m,
                  ", kmax = ",k_max,", ns = ",nstart,", mi = ",iter.max)
    if (verbose) {message("Found saved kms (",set,")...")}
    return(readRDS(fn))} else {
      
      if (verbose) {message("Calling get_cov_chol...")}
      cov_chol <- get_cov_chol(band, d, m, save, data_file, region, verbose)
      if (verbose) {message("Calling cov_chol_to_kms...")}
      kms <- cov_chol_to_kms(cov_chol, k_max, nstart, iter.max, verbose)
      
      if (save) {
        saveRDS(kms, fn)
      }
      return(kms)
    }
}

# Here d is needed to determine the label set used for attribute source
get_kms_classes_brick <- function(kms, d, save = T, 
                                  data_file = "SAR_app_data", 
                                  region = "MG", 
                                  verbose = T) {
  classes <- fitted(kms, method = "classes")
  labs <- get_pix_patched("labs", d, save, data_file, region, verbose)
  classes <- as.matrix(classes)
  classes_SAR <- SAR_matrix(classes, attr_src = labs)
  return(matrix_to_brick(classes_SAR))
}

# Plot the raster for kms output classes, and write to a raster.
# Return the kms.
write_kms_plot <- function(band, k, d, m, save = T,
                     data_file = "SAR_app_data", region = "MG",
                     k_max = 50, nstart = 20, iter.max = 100,
                     save_name = "test.tif", verbose = T,
                     col = colorspace::rainbow_hcl(k),
                     plot_only = F) {
  kms <- get_kms(band, d, m, save, 
                 data_file, region, 
                 k_max, nstart, iter.max, 
                 verbose)[[k]]
  kmsb <- get_kms_classes_brick(
    kms, d, save = T, data_file = "SAR_app_data", region = "MG", verbose = T)
  pdf(file="test.pdf",width=5,height=4)
  raster::plot(kmsb[[1]], col = col)
  dev.off()
  if (!plot_only) {raster::writeRaster(kmsb, filename = save_name, overwrite = T)}
  return(kms)
}

SARGDE_threshold_grid <- function(compare_to, d = 9, 
                                  from = 9, to = 13, by = 0.1) {
  if (2 %in% compare_to) {compare_to <- compare_to - 1}
  sargdep <- get_pix_patched(band = "sargde", d = d)
  grid <- seq(from, to, by = by); results <- list()
  for (g in grid) {
    sargde_high_g <- as.matrix(as.integer(sargdep > g))
    g <- as.character(g)
    results[[g]] <- list()
    results[[g]]$sargde_high <- sargde_high_g
    results[[g]]$conf <- caret::confusionMatrix(factor(sargde_high_g), factor(compare_to))
    results[[g]]$rand <- mclust::adjustedRandIndex(sargde_high_g, compare_to)
  }
  rands <- sapply(results, function(x){x$rand})
  threshold <- grid[which.max(rands)]
  return(list(threshold = threshold, results = results))
}

cov_chol_plotly <- function(classes,
                            band = "cc", d = 9, m = 1, 
                            sample_size = 1000,
                            subset, pairs_plot = F, col = "Set1") {
  cov_chol <- get_cov_chol(band = band, d = d, m = m)
  
  if (missing(subset)) {subset <- 1:length(classes)}
  plotlydf <- data.frame(cbind(cov_chol[subset,], classes[subset]))
  samp <- sample(1:nrow(plotlydf), sample_size)
  plotlydf <- plotlydf[samp,,drop = F]
  colnames(plotlydf) <- c("c_1","c_2","c_3","class")
  
  if (pairs_plot) {
    cls <- as.vector(unique(plotlydf[,"class"]))
    X11(); pairs(plotlydf[,c("c_1","c_2","c_3")], col = plotlydf[,"class"],
                 oma=c(3,3,3,15))
    par(xpd=TRUE)
    legend("right", legend = cls,  
           fill= cls)
  }
  plotlydf$class <- as.character(plotlydf$class)
  fig <- plot_ly(plotlydf, 
                 x = ~c_1, y = ~c_2, z = ~c_3, 
                 color = ~class,
                 colors = col)
  fig <- fig %>% add_markers()
  fig
}



###########################################################################
# FUNCTIONS FOR GETTING AND CONVERTING ------------------------------------
###########################################################################

# (1) Load saved log-cholesky vectorised autocovariances of each
#     pixel in the chosen band, or else calculate them.
get_cov_chol <- function(band, d, m, save = T,  
                         data_file = "SAR_app_data", 
                         region = "MG", verbose = T) {
  fn <- paste0(data_file,"/",toupper(region),"/",
               "/cov_chol_d",d,"_m",m,"_",band,".Rds")
  if (file.exists(fn)) {
    set <- paste0("band = ", band, ", d = ", d, ", m = ", m)
    if (verbose) {message("Found saved cov_chol (",set,")...")}
    return(readRDS(fn))} else {
    
    if (verbose) {message("Calling get_covs...")}
    covs <- get_covs(band, d, m, save, data_file, region, verbose)
    if (verbose) {message("Calling covs_to_cov_chol...")}
    cov_chol <- covs_to_cov_chol(covs, m)
    
    if (save) {
      saveRDS(cov_chol, fn)
    }
    return(cov_chol)
  }
}

# (2) Load saved autocovariances of each pixel in the chosen
#     band, or else calculate them
get_covs <- function(band, d, m, save, data_file, region, verbose) {
  fn <- paste0(data_file,"/",toupper(region),"/",
               "/covs_d",d,"_m",m,"_",band,".Rds")
  if (file.exists(fn)) {
    set <- paste0("band = ", band, ", d = ", d, ", m = ", m)
    if (verbose) {message("Found saved covs (",set,")...")}
    return(readRDS(fn))} else {
    
    if (verbose) {message("Calling get_pix_patched...")}
    pix_patched <- get_pix_patched(band, d, save, data_file, region, verbose)
    if (verbose) {message("Calling pix_patched_to_covs...")}
    covs <- pix_patched_to_covs(pix_patched, m)
    
    if (save) {
      saveRDS(covs, fn) 
    }
    return(covs)
  }
}

# (3) Load patchified pixels in the chosen
#     band, or else calculate them. If band = "labs" the labels will
#     be loaded instead. 
get_pix_patched <- function(band, d, save = T, 
                            data_file = "SAR_app_data", region = "MG", 
                            verbose = T) {
  fn <- paste0(data_file,"/",toupper(region),"/",
               band,"p_d",d,".Rds")
  if (file.exists(fn)) {
    set <- paste0("band = ", band, ", d = ", d)
    if (verbose) {message(paste0("Found saved pix_patched (",set,")..."))}
    return(readRDS(fn))} else {
    
    if (verbose) {message("Calling get_pix...")}
    pix <- get_pix(band, save, data_file, region, verbose)
    if (verbose) {message("Calling pix_to_pix_patched...")}
    pix_patched <- pix_to_pix_patched(pix, d)
    
    if (save) {
      saveRDS(pix_patched, fn) 
    }
    return(pix_patched)
  }
}

# (4) Load pixels in the chosen
#     band, or else calculate them
get_pix <- function(band, save = T, 
                    data_file = "SAR_app_data", region = "MG", 
                    verbose = T) {
  fn <- paste0(data_file,"/",toupper(region),"/",
               band,".Rds")
  if (file.exists(fn)) {
    set <- paste0("band = ", band)
    if (verbose) {message("Found saved pix (",set,")...")}
    return(readRDS(fn))} else {
    
    message("Calling load_SAR_matrix...")
    band <- if (band == "labs") {"atlas"} else {band}
    if (band %in% c("atlas","sargde")) {
      nam <- paste0(band,".tif")
      ftif <- paste0(data_file,"/",toupper(region),"/", nam)
      labsb <- raster::brick(ftif)
      labs <- rsar::brick_to_matrix(labsb, drop_na = F)
      na_ind <- is.na(labs)
      if (verbose) {message(paste0("found ", sum(na_ind), 
                    " NAs in label set, converting to 0..."))}
      labs[na_ind] <- 0
      pix <- labs
    } else {
      nam <- paste0(toupper(band),"_sub_norm.tif")
      ftif <- paste0(data_file,"/",toupper(region),"/", nam)
      pix <- load_SAR_matrix(ftif)
    }
    if (save) {
      saveRDS(pix, fn) 
    }
    return(pix)
  }
}

# Calculate patched pixel SAR_matrix from pixel SAR_matrix
pix_to_pix_patched <- function(pix, d) {
  patchify_SAR_matrix(pix, d) # Legacy function name
}

# Calculate m-lag autocovariances from patched pixel SAR_matrix
pix_patched_to_covs <- function(pix_patched, m) {
  make_mlag_autocovs(pix_patched, m) # Legacy function name
}

# Convert m-lag autocovariances to log-cholesky space
covs_to_cov_chol <- function(covs, m) {
  cov_list <- split(covs,seq(nrow(covs)))
  covmat_list <- lapply(cov_list, FUN = matrix, ncol = m+1)
  covmat_list <- lapply(covmat_list, FUN = function(x){
    if (sum(x) == 0) {diag(x) <- 1}; x})
  cov_chol <- map_to_cholesky_log2(covmat_list)
  return(cov_chol)
}

###########################################################################
# LEGACY GETTING AND CONVERTING FUNCTION NAMES ----------------------------
###########################################################################

patchify_SAR_matrix <- function(m, d) {
  img_dim <- attr(m, "brick_dim")[1:2]
  
  attribs <- NULL
  patch <- function(x, fun) {
    p <- patchify(x, n = d, m = d, img_dim = img_dim, fun = fun)
    attribs <<- attributes(p)
    return(p)
  }
  m_mean <- apply(m, MARGIN = 2, FUN = patch, fun = mean)
  m_patch <- SAR_matrix( m_mean, attr_src = m, 
                         brick_ncol = attribs$padded_brick_ncol,
                         brick_nrow = attribs$padded_brick_nrow )
  attr(m_patch, "d1_pad") <- attribs$d1_pad
  attr(m_patch, "d2_pad") <- attribs$d2_pad
  return(m_patch)
}

make_mlag_autocovs <- function(pix, m) {
  covs  <- apply(t(pix), FUN = mlag_autocov, MAR = 2, m = m)
  return(t(covs))
}

# Takes a vector v and produces the m-lag autocovariance matrix of v
mlag_autocov <- function(v, m = 2) {
  res <- matrix(v, ncol = 1)
  for (l in 1:m) {
    v <- c(v[-1],NA)
    res <- cbind(res,v)
  }
  res <- res[1:(nrow(res)-m),]
  cov(res)
}

# Function takes a vector (column of a SAR_matrix) and
# aggregates within n x m patches (provided that n
# and m are divisors of img_dim[1] and img_dim[2]
# respectively). Returns a (shorter) vector of the aggregates.
# WARNING: Don't remove NA indices before patchify.
patchify <- function(v, n, m, img_dim, fun = sum) {
  d1 <- img_dim[1]
  d2 <- img_dim[2]
  
  # Padding with zeros if needed
  nrem <- d1 %% n 
  mrem <- d2 %% m
  d1_pad <- 0
  d2_pad <- 0
  if ( nrem != 0 || mrem != 0  ) {
    v1m <- matrix(v, nrow = d1, ncol = d2, byrow = T)
    if ( nrem != 0 ) {
      v1m <- rbind( v1m, matrix(0, nrow = n - nrem, ncol = d2) )
      d1_pad <- n - nrem
      d1 <- d1 + d1_pad
    }
    if ( mrem != 0 ) {
      v1m <- cbind( v1m, matrix(0, nrow = d1, ncol = m - mrem) )
      d2_pad <- m - mrem
      d2 <- d2 + d2_pad
    }
    v <- as.vector(t(v1m))
  }
  
  # All the action happens here
  a <- reticulate::array_reshape(v, dim = c(m, d2/m, n, d1/n), order = "F")
  p <- 1
  patches <- vector(length = (d1*d2)/(n*m), mode = "numeric")
  for ( k in 1:(dim(a)[4]) ) {
    for (l in 1:(dim(a)[2]) ) {
      patches[p] <- fun(a[,l,,k])
      p <- p + 1
    }
  }
  
  attr(patches, "padded_brick_nrow") <- d1/n
  attr(patches, "padded_brick_ncol") <- d2/m
  attr(patches, "d1_pad") <- d1_pad
  attr(patches, "d2_pad") <- d2_pad
  return(patches)
}

cov_chol_to_kms <- function(cov_chol, k_max, nstart, iter.max, verbose) {
  kms <- list()
  for (k in 1:k_max) {
    if (verbose) {message(paste0("Running kmeans for k = ",k))}
    kms[[k]] <- kmeans(cov_chol, centers = k, 
                       nstart = nstart, 
                       iter.max = iter.max) 
  }
  return(kms)
}


###########################################################################
# PENALTY FUNCTIONS -------------------------------------------------------
###########################################################################

## Find the maximum adjusted Rand index.

rand_index <- function(kms,labs) {
  clust <-  fitted(kms, method = "classes")
  mclust::adjustedRandIndex(clust,labs)
}

## BIC and AIC came from here:
## https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
## A method proposed by Sherry Towers http://sherrytowers.com/2013/10/24/k-means-clustering/

kmsAIC <- function(kms) {
  v <- ncol(kms$centers)   # dimension of chol vectors
  k <- nrow(kms$centers)   # number of clusters
  ssq <- kms$tot.withinss  # total within sum of squares
  return(ssq + 2*v*k)
}

kmsBIC <- function(kms) {
  v <- ncol(kms$centers)   # dimension of chol vectors
  n <- length(kms$cluster) # number of observations
  k <- nrow(kms$centers)   # number of clusters
  ssq <- kms$tot.withinss  # total within sum of squares
  return(ssq + 2*log(n)*v*k) # higher penalty than usual
}

# Fixed logLik.kmeans from library(stackoverflow) 
logLik_kmeans <- function(object, ...) structure(
  -object$tot.withinss,
  nobs = length(object$cluster),
  df = nrow(object$centers) * ncol(object$centers),
  class = 'logLik'
)

# kmsBIC2 <- function(kms) {
#   BIC(logLik_kmeans(kms)) 
# }

# kmsBIC3 <- function(kms) {
#   v <- ncol(kms$centers)   # dimension of chol vectors
#   n <- length(kms$cluster) # number of observations
#   k <- nrow(kms$centers)   # number of clusters
#   ssq <- kms$tot.withinss  # total within sum of squares
#   return(n*log(ssq) + 10*log(n)*v*k)
# }
