---
title: "NE Method"
date: 2019-12-30
tags: ["weak", "methods"]
draft: false
---

## 1. Introduction

Instead of assuming `$L^{'}L/N \rightarrow \Sigma_{L}$` for strong factors, it is assumed that `$L^{'}L \rightarrow \Sigma_{L}$` as $N, T \rightarrow \infty$ for weak factors. The results from random matrix theory (RMT) show that, even for white noise case $\Sigma_{e} = \sigma^{2}I_{N}$, PCA estimators of the loadings and factors are inconsistent as $N, T \rightarrow \infty$.

Specifically, there exists a phase transition phenomenon in the limit: if the $k$-th largest eigenvalue of population covariance matrix `${\Sigma}_{Y}$` less than the threshold $(\sqrt{N/T}+1) \sigma^{2}$, it has little chance to detect of the $k$-th factor using PCA or MLE as $T, N \rightarrow \infty$.  Define the number of detectable factors as `$\# \{ i \leq N : \xi_{i} >  (\sqrt{N/T}+1) \sigma^{2} \}$`, where $\xi_{i}$ is the $i$-th largest eigenvalue of ${\Sigma}_{Y}$, then one goal is to estimate the number of detectable factors.

The NE method is designed to estimate the number of detectable factors in white noise for high dimensional data using the results from RMT, which studies the distribution of the sample eigenvalues. The idea of the NE method is that they use the distributional properties of the factor-free $(r=0)$ sample eigenvalues to approximate the distributional properties of the $N-r$ sample eigenvalues in $\hat{\Sigma}_{Y}$, assuming the number of factors is $r$ and $r \ll N$. The NE estimator is defined as

`$$
\hat{r}_{\mathrm{SE}} =\underset{r}{\arg \min }\left\{\frac{\beta}{4}\left[\frac{T}{N}\right]^{2} t_{r}^{2}\right\}+2(r+1) \text { for } r \in \mathbb{N}: 0 \leq r<\min (N, T),
$$`

where

`$$
t_{r} =\left[(N-r) \frac{\sum_{i=r+1}^{N} {\lambda}_{i}^{2}}{\left(\sum_{i=r+1}^{N} {\lambda}_{i}\right)^{2}}-\left(1+\frac{N}{T}\right)\right] N- \left(\frac{2}{\beta} - 1\right)\frac{N}{T}.
$$`

For the details of this method, please refer [Nadakuditi and Edelman, 2008](https://ieeexplore.ieee.org/document/4493413).

## 2. R codes

Based on this algorithm, we can write the R code for the BCV method as:

{{< highlight r >}}

NE <- function(X){

  N <- nrow(X);
  n <- ncol(X);

  sv <- svd(X)$d^2/n
  sv2 <- sv^2

  M <- min(N,n)

  tvalues <- sapply(0:M, function(k) {
    t <- ((N - k) * sum(sv2[(k + 1):M])/(sum(sv[(k + 1):M])^2)
          - (1 + N/n)) * N - N/n
    return(t)
  })

  k_hat <- which.min(n/N^2 * tvalues^2 / 4 + 2 * (0:M + 1)) - 1

  return(k_hat)
}
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

The corresponding R code of applying the NE method for this DGP is:

{{<highlight r>}}

S <- 10                   # Number of Simulation
N_set <- c(50, 100)
T_set <- c(50, 100, 200)
NT_comb <- expand.grid(N_set, T_set)  # Combination of N and T pairs

r_hat_NE <- matrix(NA, S, 4)
Factor_NE <- matrix(NA,nrow(NT_comb),6)
colnames(Factor_NE) <- c("N","T", "white", "serial", "cross", "gamma")

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
    # in the DGP:
    r_hat_NE[s,1] <-  NE(X)
    r_hat_NE[s,2] <-  NE(X1)
    r_hat_NE[s,3] <-  NE(X2)
    r_hat_NE[s,4] <-  NE(X3)
  }
    Factor_NE[i,] <- c(N,T,colMeans(r_hat_NE))
}
{{< /highlight >}}

The results of the NE method for estimating the number of strong factors with different types of error terms are:

{{<highlight r>}}
       N   T white serial cross gamma
[1,]  50  50     5   12.3  13.2   6.0
[2,] 100  50     5   17.8  18.2   6.2
[3,]  50 100     5   13.2  15.1   6.0
[4,] 100 100     5   22.9  22.4   6.0
[5,]  50 200     5    9.8  16.9   6.0
[6,] 100 200     5   22.1  25.6   6.0
{{< /highlight >}}

From the results above, we can see that the NE method is not robust when when the error terms in the factor model are high serially and cross-sectionally correlated, or have non-gaussian distribution.

## 4. Weak factor estimation

The NE method is designed to estimate the number of detectable weak factors with white noise in high dimensional data. To show this, we can apply the NE method to <a href="#DGP_weak"> the DGP </a> that we used to generate weak factors before:

{{<highlight r>}}
r_hat_NE <- matrix(NA, S, 1)
Factor_NE_weak <- matrix(NA,nrow(NT_comb),3)
colnames(Factor_NE_weak) <- c("n","T","r_hat_NE")

for(i in 1:nrow(NT_comb)){
  for(s in 1:S){

    # Data generating process
    N <- NT_comb[i,1]
    T <- NT_comb[i,2]
    X <- DGP_weak(N,T,5)$Y

    # Apply the NE method to estimate the number of strong factors
    # in the DGP:    
    r_hat_NE[s,] <-  NE(X)
  }
    Factor_NE_weak[i,] <- c(N,T,mean(r_hat_NE))
}
{{< /highlight >}}

The goal is to estimation the number of strong and useful weak factors in the DGP. The results of the NE method to estimate the number of weak factors in the DGP are:

{{<highlight r>}}
       N   T r_hat_NE
[1,]  50  50      2.4
[2,] 100  50      2.7
[3,]  50 100      2.0
[4,] 100 100      2.0
[5,]  50 200      2.0
[6,] 100 200      2.0
{{< /highlight >}}

From the results above, we can see that the NE method can estimate the number of strong and useful factors in the DGP quite well in finite samples.
