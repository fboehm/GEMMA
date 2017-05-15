#' Calculate little b hat
#' @param Xmat a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects
#' @param Vg current value of genetic covariance matrix, with dimension d by d
#' @param Ve current value of error covariance matrix, with dimension d by d
#' @param Dmat a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU
#' @param y a d by n matrix of phenotype values
#' @export

calc_b_hat <- function(Xmat, Vg, Ve, Dmat, y){
  Sigma_bb_hat <- calc_Sigmabb_hat(Xmat, Ve, Vg, Dmat)
  n_mouse <- ncol(Xmat)
  kprods <- list()
  for (l in 1:n_mouse){
    deltal <- diag(Dmat)[l]
    Vl <- deltal * Vg + Ve
    kprods[[l]] <- kronecker(Xmat[, l], solve(Vl) %*% y[, l])
  }
  # sum over list 'kprods'
  sum_kprods <- 0
  for (i in 1:length(kprods)){
    sum_kprods <- sum_kprods + kprods[[i]]
  }
  return(Sigma_bb_hat %*% sum_kprods)
}

#' Convert little b hat matrix (dimension (ad) by 1) to big B hat matrix (dimension a by d), where a is the number of allelotypes at each locus (8 in case of DO) and d is number of traits
#' @param b_hat ad by 1 matrix of stacked allelotype effect sizes
#' @param n_traits number of traits examined simultaneously. Default value is 2
#' @export

convert_b_hat <- function(b_hat, n_traits = 2){
  matrix(b_hat, nrow = n_traits, byrow = FALSE)
}


#' Calculate gl_hat in EM for multivariate LMM
#'
#' @param Xmat a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects
#' @param Vg current value of genetic covariance matrix, with dimension d by d
#' @param Ve current value of error covariance matrix, with dimension d by d
#' @param Dmat a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU
#' @param y a d by n matrix of phenotype values
#' @export

calc_gl_hat <- function(Xmat, Vg, Ve, Dmat, y, l){
  b_hat <- calc_b_hat(Xmat = Xmat, Vg = Vg, Ve = Ve, Dmat = Dmat, y = y)
  B_hat <- convert_b_hat(b_hat)
  deltal <- diag(Dmat)[l]
  Vl <- deltal * Vg + Ve
  deltal * Vg %*% solve(Vl) %*% (y[, l] - B_hat %*% Xmat[, l])

}


#' Calculate el_hat in EM for multivariate LMM
#'
#' @param yl a d by 1 matrix of phenotype values for the lth subject
#' @param B_hat a matrix of allelotype effect sizes with dimensions a by d
#' @param xl a matrix of allelotypes with dimensions d by 1 for a single subject at a single marker
#' @param gl_hat a d by 1 matrix of gl_hat values
#' @export

calc_el_hat <- function(yl, B_hat, xl, gl_hat){
  return(yl - B_hat %*% xl - gl_hat)
}


#' Calculate Sigmabb_hat for use in EM for multivariate LMM
#'
#' @param Xmat a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects
#' @param Vg current value of genetic covariance matrix, with dimension d by d
#' @param Ve current value of error covariance matrix, with dimension d by d
#' @param Dmat a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU
#' @export

calc_Sigmabb_hat <- function(Xmat, Ve, Vg, Dmat){
  n_mouse <- ncol(Xmat)
  kprods <- list()
  for (l in 1:n_mouse){
    deltal <- diag(Dmat)[l]
    Vl <- deltal * Vg + Ve
    kprods[[l]] <- kronecker(Xmat[, l] %*% t(Xmat[, l]), solve(Vl))
  }
  sum_kprods <- 0
  for (i in 1:length(kprods)){
    sum_kprods <- sum_kprods + kprods[[i]]
  }
  return(sum_kprods)
}

#' Calculate Sigmalgg_hat
#'
#' @param Xmat a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects
#' @param Vg current value of genetic covariance matrix, with dimension d by d
#' @param Ve current value of error covariance matrix, with dimension d by d
#' @param Dmat a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU
#' @param l a number from 1 to n, corresponding to the subject number
#' @export

calc_Sigmalgg_hat <- function(Xmat, Vg, Ve, Dmat, l){
  deltal <- diag(Dmat)[l]
  Vl <- deltal * Vg + Ve
  term1 <- deltal * Vg %*% solve(Vl)
  xl <- Xmat[, l]
  term2 <- kronecker(t(xl), deltal * Vg %*% solve(Vl))
  term4 <- kronecker(xl, deltal * Vg %*% solve(Vl))
  # calculate term3
  term3s <- list()
  n_mouse <- ncol(Xmat)
  for (l in 1:n_mouse){
    deltal <- diag(Dmat)[l]
    Vl <- deltal * Vg + Ve
    term3s[[l]] <- kronecker(xl %*% t(xl), deltal * Vg %*% solve(Vl))
  }
  term3s_sum <- 0
  for (i in 1:length(term3s)){
    term3s_sum <- term3s_sum + term3s[[i]]
  }
  term3 <- solve(term3s_sum)
  return(term1 + term2 %*% term3 %*% term4)
}


#' Calculate Sigmalee_hat
#'
#' @param Xmat a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects
#' @param Vg current value of genetic covariance matrix, with dimension d by d
#' @param Ve current value of error covariance matrix, with dimension d by d
#' @param Dmat a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU
#' @param l a number from 1 to n, corresponding to the subject number
#' @export

calc_Sigmalee_hat <- function(Xmat, Vg, Ve, Dmat, l){
  deltal <- diag(Dmat)[l]
  Vl <- deltal * Vg + Ve
  term1 <- deltal * Vg %*% solve(Vl)
  xl <- Xmat[, l]
  term2 <- kronecker(t(xl), solve(Vl))
  term4 <- kronecker(xl, solve(Vl))
  # calculate term3
  term3s <- list()
  n_mouse <- ncol(Xmat)
  for (l in 1:n_mouse){
    deltal <- diag(Dmat)[l]
    Vl <- deltal * Vg + Ve
    term3s[[l]] <- kronecker(xl %*% t(xl), deltal * Vg %*% solve(Vl))
  }
  term3s_sum <- 0
  for (i in 1:length(term3s)){
    term3s_sum <- term3s_sum + term3s[[i]]
  }
  term3 <- solve(term3s_sum)
  return(term1 + term2 %*% term3 %*% term4)
}



#' Update Vg
#'
#' @param Xmat a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects
#' @param Vg current value of genetic covariance matrix, with dimension d by d
#' @param Ve current value of error covariance matrix, with dimension d by d
#' @param Dmat a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU
#' @param y a d by n matrix of phenotype values
#' @export

update_Vg <- function(Xmat, Vg, Ve, Dmat, y){
  n <- ncol(y)
  summands <- list()
  b_hat <- calc_b_hat(Xmat = Xmat, Dmat = Dmat, Vg = Vg, Ve = Ve, y = y)
  B_hat <- convert_b_hat(b_hat)
  for (l in 1:n){
    deltal <- diag(Dmat)[l]
    gl_hat <- calc_gl_hat()
    Sigmalgg_hat <- calc_Sigmalgg_hat(Xmat = Xmat, Vg = Vg, Ve = Ve, Dmat = Dmat, l = l)
    summands[[l]] <- (gl_hat %*% t(gl_hat) + Sigmalgg_hat) / deltal
  }
  out <- apply(X = simplify2array(summands), MARGIN = 1:2, FUN = mean)
  return(out)
}


#' Update Ve
#'
#' @param Xmat a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects
#' @param Vg current value of genetic covariance matrix, with dimension d by d
#' @param Ve current value of error covariance matrix, with dimension d by d
#' @param Dmat a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU
#' @param y a d by n matrix of phenotype values
#' @export

update_Ve <- function(Xmat, Vg, Ve, Dmat, y){
  n <- ncol(y)
  summands <- list()
  for (l in 1:n){
    el_hat <- calc_el_hat()
    Sigmalee_hat <- calc_Sigmalee_hat()
    summands[[l]] <- el_hat %*% t(el_hat) + Sigmalee_hat
  }
  out <- apply(X = simplify2array(summands), MARGIN = 1:2, FUN = mean)
  return(out)
}


#' Calculate value of restricted likelihood for specified values of y, Vg, Ve, Xmat, and Dmat
#'
#' @param Xmat a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects
#' @param Vg current value of genetic covariance matrix, with dimension d by d
#' @param Ve current value of error covariance matrix, with dimension d by d
#' @param Dmat a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU
#' @param y a d by n matrix of phenotype values
#' @export

calc_restricted_likelihood <- function(Xmat, Vg, Ve, Dmat, y){
  d <- nrow(y) # number of phenotypes
  summands <- list()
  n_mouse <- ncol(Xmat)
  for (l in 1:n_mouse){
    term1 <- - d * log(2 * pi)
    term2 <- - log(det(Ve)) / 2
    deltal <- diag(Dmat)[l]
    term3 <- - log(deltal * Vg) / 2
    # get B_hat values
    b_hat <- calc_b_hat(Xmat = Xmat, Vg = Vg, Ve = Ve, Dmat = Dmat, y = y)
    # B_hat is the matricized version of b_hat
    B_hat <- matrix(b_hat, nrow = 2)
    gl_hat <- calc_gl_hat(Xmat = Xmat, Vg = Vg, Ve = Ve, Dmat = Dmat, y = y)
    # calculate el hat values
    el_hat <- calc_el_hat(yl = y[ , l], B_hat = B_hat, xl = Xmat[, l], gl_hat = gl_hat)
    term4 <- - t(el_hat) %*% solve(Ve) %*% el_hat / 2
    Sigmalee_hat <- calc_Sigmalee_hat(Xmat = Xmat, Vg = Vg, Ve = Ve, Dmat = Dmat, l = l)
    term5 <- - psych::tr(solve(Ve) %*% Sigmalee_hat) / 2
    term6 <- - t(gl_hat) %*% solve(deltal * Vg) %*% gl_hat / 2
    Sigmalgg_hat <- calc_Sigmalgg_hat()
    term7 <- - psych::tr(solve(deltal * Vg) %*% Sigmalgg_hat) / 2
    summands[[l]] <- term1 + term2 + term3 + term4 + term5 + term6 + term7
  }
  out <- 0
  for (i in 1:n_mouse){
    out <- out + summands[[i]]
  }
  return(out)
}

