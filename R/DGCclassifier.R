#' @import DGCPCA
NULL


#' Sample y, z pairs, y is multivariate truncated normal.
gibbs_yz <- function(init_y ,lbounds, ubounds, rotmat, lambdas, sig, burnin, num_sample){
  d <- sqrt(lambdas - sig^2)
  r <- length(d)
  p <- length(init_y)
  if (p != dim(rotmat)[1]){
    stop("The nrow of rotmat does not match the dimension of y")
  }
  if (r != dim(rotmat)[2]){
    stop("The ncol of rotmat does not match the number of lambdas")
  }
  y_sample <- matrix(nrow = num_sample, ncol = p)
  z_sample <- matrix(nrow = num_sample, ncol = r)


  y <- init_y

  # burn-in period
  for (i in 1:burnin){
    # browser()
    zmean <- (d/(d^2 + sig^2)) * (y %*% rotmat)
    zvar <- sig^2/(d^2 + sig^2)
    z <- rnorm(n = r, mean = zmean, sd = sqrt(zvar))

    ymean <- (z*d) %*% t(rotmat)

    # u <- runif(n = p, min = pnorm((lbounds - ymean)/sig), max = pnorm((ubounds - ymean)/sig))

    # y <- sig*qnorm(u) + ymean

    y <- TruncatedNormal::rtnorm(n = 1, mu = ymean, sd = sig, lb = lbounds, ub = ubounds)
  }

  # store samples to matrices
  for (i in 1:num_sample){
    # browser()
    zmean <- (d/(d^2 + sig^2)) * (y %*% rotmat)
    zvar <- sig^2/(d^2 + sig^2)
    z <- rnorm(n = r, mean = zmean, sd = sqrt(zvar))

    z_sample[i,] <- z

    ymean <- (z*d) %*% t(rotmat)

    # u <- runif(n = p, min = pnorm((lbounds - ymean)/sig), max = pnorm((ubounds - ymean)/sig))
    #
    # y <- sig*qnorm(u) + ymean

    y <- TruncatedNormal::rtnorm(n = 1, mu = ymean, sd = sig, lb = lbounds, ub = ubounds)


    y_sample[i,] <- y
  }

  return(list(z_sample = z_sample, y_sample = y_sample))


}

#' Get the lower and upper bounds for a discrete observation
get_bounds <- function(x, pmfs, cdfs){
  x <- unlist(x)
  p <- length(x)
  if (length(pmfs) != p){
    stop("The list of pmfs should be the same length of vector x")
  }

  if (length(cdfs) != p){
    stop("The list of cdfs should be the same length of vector x")
  }

  pmf_vals <- rep(0, p)
  cdf_vals <- rep(0, p)

  for (j in 1:p){
    # browser()
    pmf_vals[j] <- pmfs[[j]](x[j])
  }

  for (j in 1:p){
    cdf_vals[j] <- cdfs[[j]](x[j])
  }

  return(list(lbounds = qnorm(cdf_vals - pmf_vals),
              ubounds = qnorm(cdf_vals),
              zero_prob = any(pmf_vals == 0),
              pmf_log_sum = sum(log(pmf_vals))))

}


get_low_dim_params <- function(Rtilda, dim = r){
  p <- ncol(Rtilda)
  trace_R <- sum(diag(Rtilda))

  if (dim >= p){
    stop("dim should be less than the original dimension")
  }
  firstreigs <- RSpectra::eigs_sym(Rtilda, k = dim, which = "LA")
  lambdas <- firstreigs$values
  rotmat <- firstreigs$vectors
  sigma2 <- (trace_R - sum(lambdas))/(p - dim)
  sig <- sqrt(sigma2)
  L <- rotmat%*%diag(sqrt(lambdas - sigma2))

  return(list(lambdas = lambdas, rotmat = rotmat, sig = sig, L = L))
}


log_Phi_diff <- function(u, l){
  res <- log(pnorm(u) - pnorm(l))

  if (!is.infinite(res)){
    return(res)
  } else if (l > 0 & is.infinite(log(1 - pnorm(l)))){
    return (-log(2*pi) - 1/2*l^2 - log(l) + log1p(-exp(-1/2*(u^2 - l^2))*l/u))} else if (u < 0 & is.infinite(log(pnorm(u)))) {
      return(return (-log(2*pi) - 1/2*u^2 - log(-u) + log1p(-exp(-1/2*(l^2 - u^2))*(-u)/(-l))))
    } else {
      return(dnorm((u+l)/2, log = T) + log(u - l))
    }

}


get_ratio <- function(bounds1, bounds2, L1, L2, sig1, sig2, z){
  x1u <- unlist((bounds1$ubounds - L1 %*% z)/sig1)
  x1l <- unlist((bounds1$lbounds - L1 %*% z)/sig1)
  x2u <- unlist((bounds2$ubounds - L2 %*% z)/sig2)
  x2l <- unlist((bounds2$lbounds - L2 %*% z)/sig2)
  log_numerators <- rep(0, length(x1u))
  log_denomenators <- rep(0, length(x1l))

  for (i in 1:length(log_numerators)){
    log_numerators[i] <- log_Phi_diff(x1u[i], x1l[i])
  }

  for (i in 1:length(log_denomenators)){
    log_denomenators[i] <- log_Phi_diff(x2u[i], x2l[i])
  }

  log_ratios <- log_numerators - log_denomenators

  # browser()
  # use L'Hospital to replace NaN log_ratios
  nan_idx <- (1:length(log_ratios))[is.nan(log_ratios)]

  for (i in nan_idx){
    log_ratios[i] <- dnorm((x1u[i]+x1l[i])/2, log = T) - dnorm((x2u[i]+x2l[i])/2, log = T)
  }
  # browser()

  log_ratios_sum <- sum(log_ratios)
  return(exp(log_ratios_sum))
}


get_mean_ratio <- function(bounds1, bounds2, L1, L2, sig1, sig2, zsample2){
  # browser()
  ratios <- apply(zsample2, MARGIN = 1, FUN = function(x) get_ratio(bounds1, bounds2, L1, L2, sig1, sig2, x))

  return(mean(ratios, na.rm = T))
}


get_log_Phi_diff_sum <- function(bounds, L, sig, z){
  xu <- unlist((bounds$ubounds - L %*% z)/sig)
  xl <- unlist((bounds$lbounds - L %*% z)/sig)
  log_Phi_diffs <- mapply(log_Phi_diff, u = xu, l = xl)
  return(sum(log_Phi_diffs, na.rm = T))
}

get_log_Lik <- function(bounds, L, sig, zsample, knn){
  q <- ncol(zsample)
  log_fzx <- apply(zsample, 1, function(x) {get_log_Phi_diff_sum(bounds, L, sig, x) - q/2 * log(2*pi) - sum(x^2)/2})
  Elog_fzx <- mean(log_fzx, na.rm = T)
  entropyfzx <- FNN::entropy(X = zsample, algorithm = "kd_tree", k = knn)[knn]
  return(Elog_fzx + entropyfzx)

}





invMillRatio <- function(x){
  if (x < 3) {
    return(dnorm(x)/pnorm(x, lower.tail = F))
  } else {
    return(x)
  }
}

Qratio <- function(u, l){
  if (u < l) {
    return(Qratio(l,u))
  } else {
    if (l > 3){
      return((l/u)*exp(-1/2*(u^2 - l^2)))
    } else {
      return(pnorm(u, lower.tail = F)/pnorm(l, lower.tail = F))
    }
  }
}

log_Phi_diff_der_scale <- function(u, l){
  if (u < l) {
    return(log_Phi_diff_der_scale(l, u))
  } else {
    if ((u - l) < 1e-10) {
      return(-(u+l)/2)
    } else {
      if (l >= 0) {
        Qr <- Qratio(u, l)
        if (u == Inf) {
          invRQr <- 0
        } else {
          invRQr <- invMillRatio(u)*Qr
        }
        return( (invRQr- invMillRatio(l))/(1 - Qr))
      } else if (u <= 0){
        return(-log_Phi_diff_der_scale(-l, -u))
      } else {
        return((dnorm(u) - dnorm(l))/(pnorm(u) - pnorm(l)))
      }
    }
  }
}


neg_hfun <- function(z, L, sig, bounds){
  q <- length(z)
  p <- nrow(L)
  -(1/p)*(get_log_Phi_diff_sum(bounds, L, sig, z) - (q/2)*log(2*pi) - 1/2*sum(z^2))
}

neg_hfun_gr <- function(z, L, sig, bounds){
  q <- length(z)
  p <- nrow(L)
  gr1 <- -z
  xu <- unlist((bounds$ubounds - L %*% z)/sig)
  xl <- unlist((bounds$lbounds - L %*% z)/sig)
  scales <- mapply(log_Phi_diff_der_scale, u = xu, l = xl)/sig
  # if (any(is.nan(scales))){
  #   scales[is.nan(scales)] <- sign(mapply(log_Phi_diff_der_scale, u = xu, l = xl)[is.nan(scales)])
  # }
  gr2 <- -colSums(sweep(L, MARGIN = 1, FUN = "*", STATS = scales))

  return(-(1/p)*(gr1 + gr2))
}

neg_hfun_gr_hess_scale <- function(u, l){
  if (u < l) {
    return(neg_hfun_gr_hess_scale(l, u))
  } else {
    if ((u - l) < 1e-10) {
      return(1-((u+l)/2)^2)
    } else {
      if (l >= 0) {
        Qr <- Qratio(u, l)
        if (u == Inf) {
          uinvRQr <- 0
        } else {
          uinvRQr <- u*invMillRatio(u)*Qr
        }
        return((uinvRQr - l*invMillRatio(l))/(1 - Qr))
      } else if (u <= 0){
        return(neg_hfun_gr_hess_scale(-l, -u))
      } else {
        u <- ifelse(u == Inf, .Machine$double.xmax, u)
        l <- ifelse(l == -Inf, -.Machine$double.xmax, l)
        return((u*dnorm(u) - l*dnorm(l))/(pnorm(u) - pnorm(l)))
      }
    }
  }
}




neg_hfun_hess_det <- function(z, L, sig, bounds){
  q <- length(z)
  p <- nrow(L)
  xu <- unlist((bounds$ubounds - L %*% z)/sig)
  xl <- unlist((bounds$lbounds - L %*% z)/sig)
  scales1 <- mapply(log_Phi_diff_der_scale, u = xu, l = xl)
  scales2 <- mapply(neg_hfun_gr_hess_scale, u = xu, l = xl)
  scales <- scales1^2 + scales2
  # if (any(is.nan(scales))){
  #   scales[is.nan(scales)] <- 1
  # }
  tLSL <- crossprod(L, sweep(L, MARGIN = 1, FUN = "*", STATS = scales))
  return(-2*q*log(sig) - determinant(tLSL + diag(sig^2, q), logarithm = T)$modulus)
}

get_log_Lik_Laplace <- function(bounds, L, sig, hess = "explicit", gr = T, init_z){
  q <- ncol(L)
  p <- nrow(L)
  if (missing(init_z)) {
    init_z <- rep(0, ncol(L))
  }
  if (gr) {
    grad <- function(z) neg_hfun_gr(z, L, sig, bounds)
  } else {
    grad <- NULL
  }
  # browser()
  if (hess == "explicit") {
    optim_res <- optim(init_z, fn = function(z) neg_hfun(z, L, sig, bounds), gr = grad, method = "BFGS")

    zhat <- optim_res$par

    return(-p*optim_res$value + q/2*log(2*pi) - 1/2*neg_hfun_hess_det(zhat, L, sig, bounds) - q/2*log(p))
  } else if (hess == "approx") {
    optim_res <- optim(init_z, fn = function(z) neg_hfun(z, L, sig, bounds), gr = grad, method = "BFGS", hessian = T)

    return(-p*optim_res$value + q/2*log(2*pi) - 1/2*determinant(optim_res$hessian, logarithm = T)$modulus - q/2*log(p))
  }
}

#' Fit the Gaussian copula classification model
#'
#' @param cl The class labels of the dataset
#' @param cdfs The marginal cdf estimations used. Possible values are "empirical" (use empirical cdfs for marginal cdfs), "multinomial-dirichlet" (estimated marginal cdfs in each class shrinks to the overall marginal cdfs using multinomial-dirichlet model), or a list of estimated cdfs of class \code{stepfun}. If a list, the data are assumed to be mapped to all numerical valued data already and the support in each margin should match the support (knots) of the cdfs provided. Each item in the list is a list of marginal cdfs for class k, and the item name should match the class name in \code{cl}.
#' @param r The number of factors in the correlation matrix approximation.
#' @param ... Other parameters in \code{\link[DGCPCA]{DGCFit}}.
#' @inheritParams DGCPCA::DGCPCA
#' @export
DGCClassifier <- function(x, cl, maps = TRUE, r, cdfs = "empirical", ...){
  cl <- factor(cl)
  cl_lvs <- levels(cl)
  nk <- as.vector(table(cl)[cl_lvs])
  names(nk) <- cl_lvs
  K <- length(nk)
  fits <- list()
  if (isFALSE(maps)) {
    mapped_data <- x
  } else {
    if (isTRUE(maps)) {
      maps <- lapply(x, function(x) sort(unique(x)))
      mapped_data <- discrete_mapping(x, maps)
    } else if (is.list(maps)) {
      mapped_data <- discrete_mapping(x, maps)
    } else {
      stop("maps must be logical or a list of vectors.")
    }
  }
  if (cdfs == "empirical"){
    for (k in cl_lvs) {
      fits[[k]] <- DGCPCA::DGCFit(mapped_data[cl==k,], ...)
      fits[[k]]$low_dim_par <- get_low_dim_params(fits[[k]]$corr, dim = r)
    }
  } else if (is.list(cdfs)){
    for (k in cl_lvs) {
      fits[[k]] <- DGCPCA::DGCFit(mapped_data[cl==k,], Fhats = cdfs[[k]], ...)
      fits[[k]]$low_dim_par <- get_low_dim_params(fits[[k]]$corr, dim = r)
    }
  } else if (cdfs == "multinomial-dirichlet"){
    alpha <- list()
    supp <- list()
    for (j in 1:ncol(mapped_data)){
      # browser()
      supp[[j]] <- sort(unique(mapped_data[,j]))
      nkc <- table(mapped_data[,j], cl)
      init_a <- as.vector(table(mapped_data[,j])/nrow(mapped_data))
      init_a <- init_a * length(init_a)
      obj_f <- function(a){
        K*sum(lgamma(a)) - K*lgamma(sum(a)) + sum(lgamma(sum(a) + nk)) - sum(lgamma(a + nkc))}
      alpha[[j]] <- optim(par = init_a, fn = obj_f, method = "L-BFGS-B", lower = 1e-5)$par
      names(alpha[[j]]) <- rownames(nkc)
      # probs <- probs/sum(probs)
      # overall_cdfs[[j]] <- stepfun(x = as.numeric(rownames(nkc)), y = c(0, cumsum(probs)), right = FALSE)
    }

    for (k in cl_lvs) {
      cdfs <- list()
      # browser()
      for (j in 1:ncol(mapped_data)){
        # cats <- sort(unique(mapped_data))
        suppkj <- sort(unique(mapped_data[cl==k, j]))
        nc <- table(mapped_data[cl == k, j])[match(supp[[j]], suppkj)]
        nc[is.na(nc)] <- 0
        nc <- as.vector(nc)
        alpha_pos <- nc + alpha[[j]]
        probs <- alpha_pos/sum(alpha_pos)
        cdfs[[j]] <- stepfun(x = supp[[j]], y = c(0, cumsum(probs)[-length(probs)], 1), right = FALSE)
      }
      fits[[k]] <- DGCPCA::DGCFit(mapped_data[cl==k,], cdfs, ...)
      fits[[k]]$low_dim_par <- get_low_dim_params(fits[[k]]$corr, dim = r)
    }
  } else {
    stop("cdfs must be empirical, multinomial-dirichlet or a list")
  }

  res <- list()
  res$fits <- fits
  res$mapped_data <- mapped_data
  res$cl <- cl
  res$r<- r
  class(res) <- "DGCClassifier"
  return(res)

}

#' S3 prediction method for DGCClassifier
#'
#' @param x A DGCClassifier object
#' @param newdata The newdata to be classified. If missing, the training data will be predicted.
#' @param burin The number if burn-in period in Gibbs sampling.
#' @param num_sample The number of samples sampled in Gibbs sampling.
#' @param method method used to estimate likelihood ratio.
#' @param base.class The base class used if method = "ratio". Possible values are "max_margin", "average" or a fixed class.
#' @param knn The number of nearest neighbors to search if method = "entropy".
#' @param gr Logical. If FALSE, the gradient is approximated use finite differentiation.
#' @param hess The method to calculate the log determinant if the hessian matrix. "explicit" uses explicit formula and "approx" uses hessian approximated by finite differentiation in BFGS algorithm.
#' @return A list with a vector of predicted class and a matrix of predicted posterior probability.
#' @export
predict.DGCClassifier <- function(x, newdata, method = "ratio", base.class = "average", knn = 10, gr = T, hess = "explicit", burnin = 300, num_sample = 3000, ...){
  if (missing(newdata)) {
    newdata <- x$mapped_data
  }

  res <- list()
  res$class <- rep(NA, nrow(newdata))
  res$probs <- matrix(0, nrow = nrow(newdata), ncol = nlevels(x$cl))
  colnames(res$probs) <- levels(x$cl)

  res$mk <- rep(NA, nrow(newdata))

  cl_lvs <- levels(x$cl)
  K <- length(cl_lvs)
  pik <- table(x$cl)/length(x$cl)


  if (method == "ratio"){
    for (i in 1:nrow(newdata)){
      # browser()
      bounds <- lapply(cl_lvs, function(k) get_bounds(newdata[i,], pmfs = DGCPCA:::pmf_hat(x$fits[[k]]$cdfs), cdfs = x$fits[[k]]$cdfs))
      names(bounds) <- cl_lvs

      zero_probs <- sapply(cl_lvs, function(x) bounds[[x]]$zero_prob)
      pmf_log_sums <- sapply(cl_lvs, function(x) bounds[[x]]$pmf_log_sum)
      # browser()
      res$mk[i] <-cl_lvs[which.max(pmf_log_sums)]

      nonzero_cl <- cl_lvs[!zero_probs]

      if (length(nonzero_cl) == 0) {
        # In this case this newdata[i,] cannot be classified.
        res$probs[i, ] <- rep(0, K)
      } else if (length(nonzero_cl) == 1) {
        res$class[i] <- nonzero_cl
        res$probs[i,] <- rep(0, K)
        res$probs[i, nonzero_cl] <- 1
      } else {
        # more than 1 class has nonzero prob
        post_all <- NULL
        if (base.class == "max_margin") {
          base.class <- cl_lvs[which.max(pmf_log_sums)]
        } else if (base.class == "average") {
          base.class <- nonzero_cl
        } else if (as.character(base.class) %in% cl_lvs) {
          base.class <- as.character(base.class)
        } else {
          stop("Invalid base.class input!")
        }
        for (base_cl in base.class){
          ratios <- rep(0, K)
          names(ratios) <- cl_lvs
          ratios[base_cl] <- 1

          init_y <- qnorm(pnorm(bounds[[base_cl]]$lbounds) + runif(n = ncol(newdata)) * (pnorm(bounds[[base_cl]]$ubounds) - pnorm(bounds[[base_cl]]$lbounds)))

          yz_samples_base <- gibbs_yz(init_y = init_y, lbounds = bounds[[base_cl]]$lbounds, ubounds = bounds[[base_cl]]$ubounds, rotmat = x$fits[[base_cl]]$low_dim_par$rotmat, sig = x$fits[[base_cl]]$low_dim_par$sig, lambdas = x$fits[[base_cl]]$low_dim_par$lambdas, burnin = burnin, num_sample = num_sample)
          for (k in nonzero_cl[nonzero_cl != base_cl]){
            ratios[k] <- get_mean_ratio(bounds[[k]], bounds[[base_cl]], x$fits[[k]]$low_dim_par$L, x$fits[[base_cl]]$low_dim_par$L, x$fits[[k]]$low_dim_par$sig, x$fits[[base_cl]]$low_dim_par$sig, zsample2 = yz_samples_base$z_sample)

          }
          # browser()
          if (any(ratios == Inf)){
            ratios[ratios == Inf] <- .Machine$double.xmax
          }
          post_k <- pik*ratios/sum(pik*ratios)
          post_all <- rbind(post_all, post_k)

          # res$probs[i, base_cl] <- pik[base_cl]/sum(ratios*pik)
        }

        # mean_log_diff <- apply(log(ratios_all), 2, mean)
        # base_cl <- cl_lvs[which.max(pmf_log_sums)]
        # subractor <- min(mean_log_diff[is.finite(mean_log_diff)])
        # res$probs[i,] <- apply(post_all*as.vector(pik[nonzero_cl]), 2, sum)
        res$probs[i,] <- apply(post_all, 2, mean)
        res$probs[i,] <- res$probs[i,]/sum(res$probs[i,])
        res$class[i] <- cl_lvs[which.max(res$probs[i,])]
        # basek[i] <- base_cl
      }
      #

    }
  }
  # res$basek <- basek
  if (method == "entropy"){
    res$logLikk <- matrix(-Inf, nrow = nrow(newdata), ncol = K)
    colnames(res$logLikk) <- cl_lvs
    for (i in 1:nrow(newdata)){
      # browser()
      bounds <- lapply(cl_lvs, function(k) get_bounds(newdata[i,], pmfs = DGCPCA:::pmf_hat(x$fits[[k]]$cdfs), cdfs = x$fits[[k]]$cdfs))
      names(bounds) <- cl_lvs

      zero_probs <- sapply(cl_lvs, function(x) bounds[[x]]$zero_prob)
      pmf_log_sums <- sapply(cl_lvs, function(x) bounds[[x]]$pmf_log_sum)
      # browser()
      res$mk[i] <-cl_lvs[which.max(pmf_log_sums)]

      nonzero_cl <- cl_lvs[!zero_probs]

      if (length(nonzero_cl) == 0) {
        # In this case this newdata[i,] cannot be classified.
        res$probs[i, ] <- rep(0, K)
      } else if (length(nonzero_cl) == 1) {
        res$class[i] <- nonzero_cl
        res$probs[i,] <- rep(0, K)
        res$probs[i, nonzero_cl] <- 1
      } else {
        # more than 1 class has nonzero prob
        for (k in nonzero_cl){

          init_y <- qnorm(pnorm(bounds[[k]]$lbounds) + runif(n = ncol(newdata)) * (pnorm(bounds[[k]]$ubounds) - pnorm(bounds[[k]]$lbounds)))

          yz_samples_k <- gibbs_yz(init_y = init_y, lbounds = bounds[[k]]$lbounds, ubounds = bounds[[k]]$ubounds, rotmat = x$fits[[k]]$low_dim_par$rotmat, sig = x$fits[[k]]$low_dim_par$sig, lambdas = x$fits[[k]]$low_dim_par$lambdas, burnin = burnin, num_sample = num_sample)
          res$logLikk[i, k] <- get_log_Lik(bounds[[k]], L = x$fits[[k]]$low_dim_par$L, sig = x$fits[[k]]$low_dim_par$sig, zsample = yz_samples_k$z_sample, knn = knn)

          }
          # browser()
          for (k in nonzero_cl){
            res$probs[i, k] <- pik[k]/sum(exp(res$logLikk[i,nonzero_cl] - res$logLikk[i, k])*as.vector(pik[nonzero_cl]))
          }

          # res$probs[i, base_cl] <- pik[base_cl]/sum(ratios*pik)
        }

        # mean_log_diff <- apply(log(ratios_all), 2, mean)
        # base_cl <- cl_lvs[which.max(pmf_log_sums)]
        # subractor <- min(mean_log_diff[is.finite(mean_log_diff)])
        # res$probs[i,] <- apply(post_all*as.vector(pik[nonzero_cl]), 2, sum)
        # res$probs[i,] <- apply(post_all, 2, mean)
        # res$probs[i,] <- res$probs[i,]/sum(res$probs[i,])
        res$class[i] <- cl_lvs[which.max(res$probs[i,])]
        # basek[i] <- base_cl
      }
      #

  }

  if (method == "laplace"){
    res$logLikk <- matrix(-Inf, nrow = nrow(newdata), ncol = K)
    colnames(res$logLikk) <- cl_lvs
    for (i in 1:nrow(newdata)){
      # browser()
      bounds <- lapply(cl_lvs, function(k) get_bounds(newdata[i,], pmfs = DGCPCA:::pmf_hat(x$fits[[k]]$cdfs), cdfs = x$fits[[k]]$cdfs))
      names(bounds) <- cl_lvs

      zero_probs <- sapply(cl_lvs, function(x) bounds[[x]]$zero_prob)
      pmf_log_sums <- sapply(cl_lvs, function(x) bounds[[x]]$pmf_log_sum)
      # browser()
      res$mk[i] <-cl_lvs[which.max(pmf_log_sums)]

      nonzero_cl <- cl_lvs[!zero_probs]

      if (length(nonzero_cl) == 0) {
        # In this case this newdata[i,] cannot be classified.
        res$probs[i, ] <- rep(0, K)
      } else if (length(nonzero_cl) == 1) {
        res$class[i] <- nonzero_cl
        res$probs[i,] <- rep(0, K)
        res$probs[i, nonzero_cl] <- 1
      } else {
        # more than 1 class has nonzero prob
        for (k in nonzero_cl){

          init_y <- qnorm(pnorm(bounds[[k]]$lbounds) + runif(n = ncol(newdata)) * (pnorm(bounds[[k]]$ubounds) - pnorm(bounds[[k]]$lbounds)))

          yz_samples_k <- gibbs_yz(init_y = init_y, lbounds = bounds[[k]]$lbounds, ubounds = bounds[[k]]$ubounds, rotmat = x$fits[[k]]$low_dim_par$rotmat, sig = x$fits[[k]]$low_dim_par$sig, lambdas = x$fits[[k]]$low_dim_par$lambdas, burnin = 0, num_sample = 1)
          res$logLikk[i, k] <- get_log_Lik_Laplace(bounds = bounds[[k]], L = x$fits[[k]]$low_dim_par$L, sig = x$fits[[k]]$low_dim_par$sig, hess = hess, init_z = yz_samples_k$z_sample[1,])

        }
        # browser()
        for (k in nonzero_cl){
          res$probs[i, k] <- pik[k]/sum(exp(res$logLikk[i,nonzero_cl] - res$logLikk[i, k])*as.vector(pik[nonzero_cl]))
        }

        # res$probs[i, base_cl] <- pik[base_cl]/sum(ratios*pik)
      }


      res$class[i] <- cl_lvs[which.max(res$probs[i,])]
    }
    #

  }

  return(res)
}
