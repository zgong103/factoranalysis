---
title: "BCV Method"
date: 2019-12-30
tags: ["weak", "methods"]
draft: false
---

## 1. Introduction

Instead of estimating the number of detectable factors, one may prefer estimating the number of useful factors (including strong and useful weak factors). The number of useful factors recover an underlying signal matrix $X = LF$ in the factor model more precisely than using the true number of factors or detectable factors. The BCV method is designed to estimate the number of useful factors in heteroscedastic noise for high dimensional data based on bi-cross-validation, using randomly held-out submatrices of the data matrix.

The algorithm for recovering the signal matrix $X$ has two steps. First, they devise early stopping alternation (ESA) method to estimate $X$ given the number of factors $r^{* }$. Second, they propose bi-cross-validation (BCV) method to estimate the number of factors $r^{* }$ based on the ESA method. The idea of BCV method is that, for each candidate $r$,  we first use the three held-in blocks to estimate the held-out block $X_{00}$ (corresponding to $Y_{00}$ in the factor model) and then select the optimal $r^*$ by minimizing the BCV estimated prediction error. For the details of this method, please refer [Owen and Wang, 2015](https://arxiv.org/abs/1503.03515).



## 2. R codes

Based on the algorithm above, we can write the R code for the BCV method as:

{{< highlight r >}}

EsaBcv <- function (Y, X = NULL, r.limit = 20, niter = 3, nRepeat = 12,
                    only.r = F, svd.method = "fast", center = F){

    Y <- as.matrix(Y);
    p <- ncol(Y);
    n <- nrow(Y);
    X.ori <- X;

    # center is a logical term: whether to add an intercept term in the model.
    if (center)         
    X <- cbind(rep(1, n), X);
    qr.X <- NULL;
    Y.ori <- Y;

    if(!is.null(X)) {
		X <- as.matrix(X);
		k <- ncol(X);
		if (k >= n)
			stop("Too many predictors!! The number of predictors
			     is expected to be much less than the sample size!")
		qr.X <- qr(X);
		Y <- qr.qty(qr.X, Y)[-(1:k), ];
		n <- n-k;
    }

    ## decide the held-in size
    gamma <- p/n;
    bar.gamma <- ((sqrt(gamma) + sqrt(1/gamma))/2)^2;
    sqrt.rho <- sqrt(2)/(sqrt(bar.gamma) + sqrt(bar.gamma + 3));
    held.in.size.small <- min(round(sqrt.rho*sqrt(p * n)), p - 1, n - 1);
    held.in.size.large <- round(sqrt.rho^2 * p * n/held.in.size.small);
    if (n < p) {
        n1 <- held.in.size.small;
        p1 <- held.in.size.large;
    } else {
        n1 <- held.in.size.large;
        p1 <- held.in.size.small;
    }
    max.r <- min(n1, p1);
    if(max.r > r.limit)
        max.r <- r.limit;
    if (is.null(nRepeat)) {
       nRepeat <- max(round(p/p1), round(n/n1), round(p/(p-p1)),
                      round(n/(n-n1)));
    }

    result.list <- NULL;

   # Partition of the data matrix and choose r with the smallest prediction error
    for (i in 1:nRepeat) {
       	Y.resample <- Y[sample(1:n, n), sample(1:p, p)];
        result.list <- rbind(result.list,sapply(0:max.r, function(r)
          BcvGivenRank(Y.resample, r, niter,n1,p1,svd.method)));
    }
    not.na.result <- which(!is.na(colMeans(result.list)));
    max.r <- sum(not.na.result == 1:length(not.na.result)) - 1;
    result.list <- result.list[, 1:(max.r + 1)];
	  colnames(result.list) <- 0:max.r;
    best.r <- as.numeric(which.min(colMeans(result.list)) - 1);

    if (only.r)
    return(best.r)

    est <- ESA(Y.ori, best.r, X.ori, center, niter, svd.method);
    result <- list(best.r = best.r, estSigma = est$estSigma, estU = est$estU,     estD = est$estD, estV = est$estV, beta = est$beta,
		estS = est$estS, mu = est$mu, max.r = max.r);
	  class(result) <- "esabcv"
    return(result);
}
{{< /highlight >}}

In the function above, we return a vector of BCV MSE from the four folds for a given rank:
{{<highlight r>}}
BcvGivenRank <- function (Y, r, niter, n1, p1, svd.method) {
    p <- ncol(Y);
    n <- nrow(Y);
    hp <- p - p1;
    hn <- n - n1;
    bcv.result <- BcvPerFoldGivenRank(Y[1:hn, 1:hp, drop = F],
                                      Y[1:hn,-(1:hp), drop = F],
                                      Y[-(1:hn),1:hp, drop = F],
                                      Y[-(1:hn),-(1:hp), drop = F], r, niter,                                       svd.method);
    return(bcv.result);
}
{{< /highlight >}}

We calculate BCV MSE for one fold:

{{<highlight r>}}
## inputs:
# Y00, Y01, Y10, Y11: the four folds of data and Y00 is the held-out fold
# r: given number of factor to try
# niter: number of iteration steps for ESA.
## outputs:
# the average entrywise estimation error of the held-out block.
# if the estimate of Sigma is unreasonable, return NA

BcvPerFoldGivenRank <- function(Y00, Y01, Y10, Y11, r, niter, svd.method,
                                  tol.Sigma.scale = 6){
    if (r ==0)
        return(mean(Y00^2));
    result <- try(ESA(Y11, r, niter = niter, svd.method = svd.method));
    if (class(result)== "try-error") {
        save(Y11, r, file = "problematic.RData")
        print("Encounter problematic data!")
    }
    Sigma1 <- result$estSigma;

    if ((mean(-log10(Sigma1)) > tol.Sigma.scale - log10(max(diag(Sigma1)))) |
        (max(diag(Sigma1)) < .Machine$double.eps)) {
        return(NA);
    }

    Y11.tilde <- t(t(Y11)/ sqrt(Sigma1));
    Y01.tilde <- t(t(Y01) / sqrt(Sigma1));
    pseudo.inv.Y11.tilde <- PseudoInv(Y11.tilde, r, svd.method);
    Y00.est <- Y01.tilde %*% pseudo.inv.Y11.tilde %*% Y10;
    held.out.res <- Y00 - Y00.est;
    est.error <-  mean((held.out.res)^2);

    return (est.error);
}
{{< /highlight >}}

We calculate the Moore-Penrose pseudo inverse of matrix Y up to a given
rank, using moore-penrose pseudo inverse of a matrix as:

{{<highlight r>}}
# Y the given matrix
# k the given rank
# A pseudo-inverse matrix of rank \code{k}
PseudoInv <- function(Y, k, svd.method = "fast") {
      svd.trunc <- ComputeSVD(Y, svd.method, k);
      tol <- sqrt(.Machine$double.eps);
      pseudo.inv.d <- sapply(svd.trunc$d, function(x)
                             ifelse(x > svd.trunc$d[1] * tol, x^(-1), 0));
      pseudo.inv.Y <- svd.trunc$v%*%diag(pseudo.inv.d[1:k],k,k)%*%
                                                  t(svd.trunc$u);
      return(pseudo.inv.Y)
}
{{< /highlight >}}

We estimate the latent factor matrix and noise variance using early stopping alternation (ESA) given the number of factors as:

{{<highlight r>}}
ESA <- function(Y, r, X = NULL, center = F, niter = 3, svd.method = "fast"){

    Y <- as.matrix(Y);
    p <- ncol(Y);
    n <- nrow(Y);
    Y.ori <- Y;
    if (center)
    	X <- cbind(rep(1, n), X);
    qr.X <- NULL;
    if(!is.null(X)) {
		X <- as.matrix(X);
		k <- ncol(X);
		if (k >= n)
			stop("Too many predictors!! The number of predictors
			     is expected to be much less than the sample size!")
		qr.X <- qr(X);
		beta <- qr.coef(qr.X, Y);
		if (center) {
			mu <- beta[1, ];
			if (k == 1) {
				beta1 <- NULL;
			} else
				beta1 <- beta[-1, ];
		} else {
			mu <- NULL;
			beta1 <- beta;
		}
		Y <- qr.qty(qr.X, Y)[-(1:k), ];
		n <- n-k;
    } else {
        beta1 <- NULL;
        mu <- NULL;
    }

    # initializing Sigma
    Sigma <- apply(Y, 2, var);
    if (r == 0)
        return(list(estSigma = Sigma, estU = NULL, estD = NULL,
                    estV = NULL, estS = NULL, beta = beta1, mu = mu));
    if (r >= min(p, n)) {
        svd.Y <- ComputeSVD(Y, svd.method = svd.method);
        if (!is.null(X)) {
        	estU <- qr.qy(qr.X, rbind(matrix(rep(0, k * r), nrow = k),
        						      svd.Y$u))
        } else
        	estU <- svd.Y$u
        return(list(estSigma = rep(0,p), estU = estU,
                    estD = svd.Y$d / sqrt(n),
                    estV = svd.Y$v, estS = Y.ori, beta = beta1, mu = mu));   
    }
    iter <- 0;
    while (1) {
        iter = iter + 1;
        if( iter > niter ){
            break;
        }
        Y.tilde <- sapply(1:p, function(i)
                            Y[, i] / sqrt(Sigma[i]));
        svd.Ytilde <- ComputeSVD(Y.tilde, rank = r);
        U.tilde <- svd.Ytilde$u;
        V.tilde <- svd.Ytilde$v;
        estU <- U.tilde;
        estD <- diag(svd.Ytilde$d[1:r], r, r)
        estV <- sqrt(Sigma) * V.tilde;     
        res <- Y - estU %*% estD %*% t(estV);
        Sigma <- apply(res, 2, function(x) sum(x^2)/length(x));
    }
    if (!is.null(X))
		estU <- rbind(matrix(rep(0, k * r), nrow = k), estU);
    estS <- estU %*% estD %*% t(estV);
    svd.Y <- ComputeSVD(estS, svd.method, r);
	if (!is.null(X)) {
		estU <- qr.qy(qr.X, svd.Y$u);
		estS <- estS + X %*% beta;
	} else
		estU <- svd.Y$u;
    estD <- svd.Y$d[1:r]/sqrt(n);
    return(list(estSigma = Sigma, estU = estU, estD = estD,
                estV = svd.Y$v, estS = estS, beta = beta1, mu = mu));    
}
{{< /highlight >}}

A wrapper for computing SVD:

{{<highlight r>}}
ComputeSVD <- function(Y, svd.method = "fast", rank = NULL, kmax.r = 10) {
    if (is.null(rank)) {
        rank <- min(dim(Y));
    }
	if(svd.method == "propack") {
		svd.Y <- try(suppressWarnings(propack.svd(Y, neig = rank)));
		if (class(svd.Y) == "try-error" & (rank < min(dim(Y)))) {
			svd.Y <- suppressWarnings(propack.svd(Y, neig = rank + 1));
			if (!length(svd.Y$d) == (rank + 1)) {
				svd.Y <- propack.svd(Y, neig = rank + 1,
									 opts = list(kmax = kmax.r * (rank + 1)));
			}
			svd.Y$u <- svd.Y$u[, 1:rank];
			svd.Y$v <- svd.Y$v[, 1:rank];
			svd.Y$d <- svd.Y$d[1:rank];
		} else if ((class(svd.Y) == "try-error") | (length(svd.Y$d) != rank)) {
			warning("PROPACK fails, used fast.svd to compute SVD instead!!");
			svd.Y <- ComputeSVD(Y, "fast", rank);
		}
	} else if (svd.method == "fast") {
        tol <- max(dim(Y)) * .Machine$double.eps;

        svd.Y <- fast.svd(Y, tol);

        if (rank < min(dim(Y))) {
            svd.Y$u <- matrix(svd.Y$u[, 1:rank], ncol = rank);
            svd.Y$v <- matrix(svd.Y$v[, 1:rank], ncol = rank);
            svd.Y$d <- svd.Y$d[1:rank];
        }
    } else    {
        svd.Y <- svd(Y, nu = rank, nv = rank);
    }
    return(svd.Y)
}
{{< /highlight >}}

The package we need to `library()` before applying the BCV method is:

{{<highlight r>}}
library(corpcor)
{{< /highlight >}}




## 3. Strong factor estimation and robustness check

Simulation design: we use the following DGP </a> for generating strong factors:

\begin{equation}
    Y_{it} = \sum_{j=1}^{r}\lambda_{ij} F_{tj}  + e_{it}, \quad \text{where}\\
    \lambda_{ij}, F_{tj} \stackrel{\text { iid }}{\sim} \mathcal{N}(0,1),\\
    \end{equation}
$$\text{with} \quad e_{i t} = \rho_{1} e_{i t-1} + (1-\rho_{1}^2)^{1/2} \xi_{it},$$
$$\text{and} \quad \xi_{i t} = \rho_{2} \xi_{i-1, t} + (1-\rho_{2}^2)^{1/2} \epsilon_{it}, \quad \epsilon_{it} \stackrel{\text { iid }}{\sim} \mathcal{N}(0,1).$$

We let $r=5$,  and consider the three cases for $e_{it}$ below:

  * **Case I:**   high serial correlation only, $\rho_{1} = 0.9$ and $\rho_{2} = 0$;

  * **Case II:**  high cross-sectional correlation only, $\rho_{1} = 0$ and $\rho_{2} = 0.8$;

  * **Case III:** non-guassion distributions only, $\rho_{1} = \rho_{2} = 0$. For example, we can chosse gamma distribution for $e_{it}$ with mean zero and variance 0.5.

The corresponding R code of applying the BCV method for this DGP is:

{{<highlight r>}}

S <- 10                   # Number of Simulation
N_set <- c(50, 100)
T_set <- c(50, 100, 200)
NT_comb <- expand.grid(N_set, T_set)  # Combination of N and T pairs

r_hat_BCV <- matrix(NA, S, 4)
Factor_BCV <- matrix(NA,nrow(NT_comb),6)
colnames(Factor_BCV) <- c("N","T", "white", "serial", "cross", "gamma")

for(i in 1:nrow(NT_comb)){
  for(s in 1:S){

    # Data generating process
    N <- NT_comb[i,1]
    T <- NT_comb[i,2]

    F_0 <-  matrix(rnorm(T*r), T, r)  # generating factor matrix
    L_0 <-  matrix(rnorm(N*r), N, r)  # generating loading matrix

    # DGP for generating strong factors in white noise
    e_0 <-  matrix(rnorm(N*T), N, T)  
    X <- L_0%*%t(F_0) + e_0   

    # DGP for generating strong factors in high serially correlated noise
    rho1 <- 0.8
    e_1 <- matrix(NA, N, T)
    e_1[, 1] <- rnorm(N, 0, 1)
    for (t in 1:(T-1)) {
      e_1[, t+1] <- e_1[, t]*rho1 + sqrt(1 - rho1^2)*rnorm(N, 0, 1)
    }
    X1 <- L_0%*%t(F_0) + e_1

    # DGP for generating strong factors in high cross-sectionally
    # correlated noise
    rho2 <- 0.8
    e_2 <- matrix(NA, N, T)
    e_2[1,] <- rnorm(T, 0, 1)
    for (n in 1:(N-1)) {
      e_2[n+1 , ] <- e_2[n, ]*rho2 + sqrt(1 - rho2^2)*rnorm(T, 0, 1)
    }
    X2 <- L_0%*%t(F_0) + e_2

    # DGP for generating strong factors in noise with gamma distribution
    e_3 <-  matrix(rgamma(N * T, 0.25, scale = 4), nrow = N)
    X3 <- L_0%*%t(F_0) + e_3

    # Apply the NE method to estimate the number of strong factors
    # in the DGP
    r_hat_BCV[s,1] <-  EsaBcv(X, only.r = TRUE)
    r_hat_BCV[s,2] <-  EsaBcv(X1, only.r = TRUE)
    r_hat_BCV[s,3] <-  EsaBcv(X2, only.r = TRUE)
    r_hat_BCV[s,4] <-  EsaBcv(X3, only.r = TRUE)

  }
    Factor_BCV[i,] <- c(N,T,colMeans(r_hat_BCV))
}
{{< /highlight >}}

The results of the BCV method for estimating the number of strong factors with different types of error terms are:

{{<highlight r>}}
       N   T white serial cross gamma
[1,]  50  50     5   10.0  10.6   4.9
[2,] 100  50     5   14.3  11.1   6.5
[3,]  50 100     5   11.4  13.5   6.3
[4,] 100 100     5   17.9  16.6   6.0
[5,]  50 200     5    8.7  15.5   6.1
[6,] 100 200     5   19.1  19.8   6.1
{{< /highlight >}}

From the results above, we can see that the BCV method is not robust when the error terms in the factor model are high serially and cross-sectionally correlated, or have non-gaussian distribution.

## 4. Weak factor estimation

The BCV method is designed to estimate the number of strong and useful weak factors with heteroscedastic noise in high dimensional data. To show this, we can apply the BCV method to <a href="#DGP_weak"> the DGP</a> that we used to generate weak factors before:

{{<highlight r>}}
r_hat_BCV <- matrix(NA, S, 1)
Factor_BCV_weak <- matrix(NA,nrow(NT_comb),3)
colnames(Factor_BCV_weak) <- c("n","T","r_hat_BCV")

for(i in 1:nrow(NT_comb)){
  for(s in 1:S){

# Data generating process
    N <- NT_comb[i,1]
    T <- NT_comb[i,2]
    X <- DGP_weak(N,T,5)$Y

# Apply the NE method to estimate the number of strong factors in the DGP:    
    r_hat_BCV[s,] <-  EsaBcv(X, only.r = TRUE)
  }
    Factor_BCV_weak[i,] <- c(N,T,mean(r_hat_BCV))
}
{{< /highlight >}}

The goal is to estimation the number of strong and useful weak factors in the DGP. The results of the NE method to estimate the number of weak factors in the DGP are:

{{<highlight r>}}
       N   T r_hat_BCV
[1,]  50  50       2.1
[2,] 100  50       2.8
[3,]  50 100       2.6
[4,] 100 100       3.0
[5,]  50 200       2.8
[6,] 100 200       3.0
{{< /highlight >}}

From the results above, we can see that the BCV method can estimate the number of strong and useful factors in the DGP very well in finite samples.
