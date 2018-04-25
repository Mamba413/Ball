#' @title Distance Matrix Computation for Non-Hilbert Data
#' 
#' @description This function computes and returns the numeric distance matrix 
#' computed by using the specified distance measure to compute 
#' the distances between the rows of a data matrix.
#' 
#' @param x a numeric matrix, data frame or numeric array of dimension \eqn{k \times m \times n} 
#' containing \eqn{n} samples in \eqn{k \times m} dimension.
#' @param method the distance measure to be used. This must be one of "geodesic", "compositional", "riemann". 
#' Any unambiguous substring can be given.
#'
#' @details Available distance measures are geodesic, compositional and riemann.
#' Denoting any two sample in the dataset as \eqn{x} and \eqn{y}, 
#' we give the definition of distance measures as follows.
#' 
#' geodesic:
#' 
#' The shortest route between two points on the Earth's surface, namely, a segment of a great circle.
#' \deqn{acos(x^{T}y), \|x\|_{2} = \|y\|_{2} = 1}
#' 
#' compositional:
#' 
#' First, we apply scale transformation to it, i.e., \eqn{(x_{i1}/t, ..., x_{ip}/t_{i}), t_{i} = \sum_{d=1}^{p}{x_{d}}} 
#' . Then, apply the square root transformation to data and calculate the geodesic distance between samples.
#' 
#' riemann:
#' 
#' \eqn{k \times m \times n} array where \eqn{k} = number of landmarks, \eqn{m} = number of dimensions and \eqn{n} = sample size. Detail about
#' riemannian shape distance was given in Kendall, D. G. (1984).
#' 
#' @return \eqn{n \times n} numeric distance matrix
#' @references Kendall, D. G. (1984). Shape manifolds, Procrustean metrics and complex projective spaces, Bulletin of the London Mathematical Society, 16, 81-121.
#' @export 
#' @examples
#' data('bdvmf')
#' Dmat <- nhdist(bdvmf[['x']], method = "geodesic")
#' 
#' data("ArcticLake")
#' Dmat <- nhdist(ArcticLake[['x']], method = "compositional")
#' 
#' data("macaques")
#' Dmat <- nhdist(macaques[["x"]], method = "riemann")
#' 
#' # unambiguous substring also available:
#' Dmat <- nhdist(macaques[["x"]], method = "rie")
#' 
nhdist <- function(x, method = 'geodesic') {
  METHODS <- c("geodesic", "compositional", "riemann", "bhattacharyya", "angular")
  methodIndex <- pmatch(method, METHODS)
  if (is.na(method)) 
    stop("invalid distance method")
  if (method == -1) 
    stop("ambiguous distance method")
  method <- METHODS[methodIndex]
  #
  if(method %in% c('geodesic')) {
    return(distsurface(x))
  } else if(method %in% c('compositional')) {
    return(dist_Bhattacharyya(x))
  } else if(method == "bhattacharyya") {
    return(dist_Bhattacharyya(x))
  } else if(method == "angular") {
    return(dist_Angular(x))
  } else if(method %in% c('riemann')) {
    return(distrieman(x))
  }
}

#' Distance for compositional data (Bhattacharyya distance)
#'
#' @param x Matrix object.
#'
#' @return Distance matrix
#' @noRd
#' @examples
#' data("ArcticLake")
#' head(ArcticLake)
#' distcompositional(ArcticLake[, 1:3])
dist_bhattacharyya <- function(x) {
  xRowSum <- rowSums(x)
  x <- apply(x, 2, function(z) {
    z/xRowSum
  })
  x <- sqrt(x)
  #
  distsurface(x)
}


#' Distance for compositional data (Angular distance)
#'
#' @param x Matrix object.
#'
#' @return Distance matrix
#' @noRd
#' @examples
#' data("ArcticLake")
#' head(ArcticLake)
#' distcompositional(ArcticLake[, 1:3])
dist_angular <- function(x) {
  x <- x^2
  xRowSum <- rowSums(x)
  x <- apply(x, 2, function(z) {
    z/xRowSum
  })
  x <- sqrt(x)
  #
  distsurface(x)
}



#' @title Geodesic Distance in Unit Ball
#'
#' @param x Matrix object. 
#'
#' @return Distance matrix
#' @noRd
#' @examples
#' data("bdvmf")
#' Dmat <- distsurface(bdvmf[['x']])
#' 
distsurface <- function(x) {
  x <- as.matrix(x)
  Dmat <- x %*% t(x)
  diag(Dmat) <- 1
  acos(Dmat)
}


#' @title Riemannian shape distance
#'
#' @param x \eqn{ k \times m \times n } array.
#' @description Calculates the Riemannian shape distance rho between two configurations
#' @return Distance matrix
#' @noRd
#' @examples
#' data("macaques")
#' Dmat <- distrieman(macaques[["x"]])
#' 
distrieman <- function(x) {
  n <- dim(x)[3]
  sapply(1:n, function(i) {
    sapply(1:n, function(j) {
      riemdist(x[,,i], x[,,j])
    })
  })
}

realtocomplex<-function(x)
{
  #input k × 2 matrix - return complex k-vector 
  k <- nrow(x)
  zstar <- x[, 1] + (1i) * x[, 2]
  zstar
}

defh<-function(nrow)
{
  #Defines and returns an nrow × (nrow+1) Helmert sub-matrix
  k <- nrow
  h <- matrix(0, k, k + 1)
  j <- 1
  while(j <= k) {
    jj <- 1
    while(jj <= j) {
      h[j, jj] <- -1/sqrt(j * (j + 1))
      jj <- jj + 1
    }
    h[j, j + 1] <- j/sqrt(j * (j + 1))
    j <- j + 1
  }
  h
}

st<-function(zstar)
{
  #input complex matrix
  #output transpose of the complex conjugate 
  st <- t(Conj(zstar))
  st
}


centroid.size<-function(x)
{
  #returns the centroid size of a configuration (or configurations)
  #input: k × m matrix/or a complex k-vector
  # or input a real k × m × n array to get a vector of sizes for a sample
  if ((is.vector(x)==FALSE) && is.complex(x)){
    k <- nrow(x)
    n <- ncol(x)
    tem <- array(0,c(k,2,n))
    tem[ ,1, ] <- Re(x)
    tem[ ,2, ] <- Im(x)
    x <- tem
  }
  {
    if (length(dim(x))==3){
      n <- dim(x)[3]
      sz <- rep(0,times=n)
      k <- dim(x)[1]
      h <- defh(k - 1)
      for (i in 1:n){
        xh <- h %*% x[ , ,i]
        sz[i] <- sqrt(sum(diag(t(xh) %*% xh)))      
      }
      sz
    } 
    else 
      {
      if (is.vector(x) && is.complex(x)) {
        x <- cbind(Re(x), Im(x))
      }
      k <- nrow(x)
      h <- defh(k - 1)
      xh <- h %*% x
      size <- sqrt(sum(diag(t(xh) %*% xh)))
      size
    }
  }        
}

preshape<-function(x)
{
  #input k × m matrix / complex k-vector
  #output k-1 × m matrix / k-1 × 1 complex matrix
  if(is.complex(x)) {
    k <- nrow(as.matrix(x))
    h <- defh(k - 1)
    zstar <- x
    ztem <- h %*% zstar
    size <- sqrt(diag(Re(st(ztem) %*% ztem)))
    if(is.vector(zstar))
      z <- ztem/size
    if(is.matrix(zstar))
      z <- ztem %*% diag(1/size)
  }
  else {
    if(length(dim(x)) == 3) {
      k <- dim(x)[1]
      h <- defh(k - 1)
      n <- dim(x)[3]
      m <- dim(x)[2]
      z <- array(0, c(k - 1, m, n))
      for(i in 1:n) {
        z[,  , i] <- h %*% x[,  , i]
        size <- centroid.size(x[,  , i])
        z[,  , i] <- z[,  , i]/size
      }
    }
    else {
      k <- nrow(as.matrix(x))
      h <- defh(k - 1)
      ztem <- h %*% x
      size <- centroid.size(x)
      z <- ztem/size
    }
  }
  z
}


riemdist<-function(x, y, reflect=FALSE)
{
  #input two k × m matrices x, y or complex k-vectors
  #output Riemannian distance rho between them 
  
  if (sum((x-y)**2)==0){
    riem <- 0
  }
  if (sum((x-y)**2)!=0){
    
    if (reflect==FALSE) { 
      if(ncol(as.matrix(x)) < 3) {
        if (is.complex(x)==FALSE){x<-realtocomplex(x)}
        if (is.complex(y)==FALSE){y<-realtocomplex(y)} 
        #riem <- c(acos(Mod(st(preshape(x)) %*% preshape(y))))
        riem<-c(acos(min(1,(Mod(st(preshape(x)) %*% preshape(y))))))
      }
      else {
        m <- ncol(x)
        z <- preshape(x)
        w <- preshape(y)
        Q <- t(z) %*% w %*% t(w) %*% z
        ev <- eigen(t(z) %*% w)$values
        check <- 1
        for(i in 1:m) {
          check <- check * ev[i]
        }
        ev <- sqrt(abs(eigen(Q, symmetric = TRUE)$values))
        if(Re(check) < 0)
          ev[m] <-  - ev[m]
        riem <- acos(min(sum(ev),1))
      }
    }
    if (reflect==TRUE){
      
      m <- ncol(x)
      z <- preshape(x)
      w <- preshape(y)
      Q <- t(z) %*% w %*% t(w) %*% z
      
      ev <- sqrt(abs(eigen(Q, symmetric = TRUE)$values))
      
      riem <- acos(min(sum(ev),1))
    }
    
  }
  riem
}