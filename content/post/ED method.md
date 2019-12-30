---
title: "ED Method"
date: 2019-12-29
tags: ["strong", "methods"]
draft: false
---

## 1. Introduction

The ED estimator is defined as

`$$\hat{r}_{\mathrm{ED}} = \max \left\{r \leq r_{max}: \lambda_{r} - \lambda_{r+1} \geq \delta\right\},$$`

where $\delta$ is some fixed number, $\lambda_{i}$ is the $i$-th largest eigenvalue of  $\hat{\Sigma}_{Y}$. This method estimates the number of factors by exploiting the structure of idiosyncratic terms using the results from RMT. It explicitly allows serial and cross-sectional correlation in the error terms in the assumption.

An advantage of this method comparing with the IC method is that the consistency of the ED estimator can allow for much weaker strength of the factors: instead of growing in the order of $O(N)$, the smallest eigenvalue of $L'L$ are just required to diverge in probability as $N \rightarrow \infty$. The algorithm for choosing $\delta$ and estimating the number of factors is:

1. Compute eigenvalues $\lambda_{1}, \ldots, \lambda_{n}$ of the sample covariance matrix $X X^{\prime} / T$. Set $j=r_{\text {max }}+1 .$
2. Compute $\hat{\beta}$, the slope coefficient in the OLS regression
of $\lambda_{j}, \ldots, \lambda_{j+4},$ on the constant and $(j-1)^{2 / 3}, \ldots,(j+3)^{2/3}.$ Set $\delta=2|\hat{\beta}|$.
3. Compute $\hat{r}(\delta) = \max \{\lambda_{i}-\lambda_{i+1} \geq \delta \}$ or if $\lambda_{i}-\lambda_{i+1} < \delta$ for all $i \leq r_{\text {max }}^{n},$ set $\hat{r}(\delta)=0$.
4. Set $j=\hat{r}(\delta)+1 .$ Repeat steps 2 and 3 until conver-
gence.

For the details of this method, please refer [Onatski, 2010](https://www.mitpressjournals.org/doi/abs/10.1162/REST_a_00043?journalCode=rest).

## 2. R codes

Based on the algorithm above, we can write the R code for the ED method as:

{{< highlight r >}}

EigenDiff <- function(Y, rmax = 20, niter = 10) {

  N <- nrow(Y)
  T <- ncol(Y)

  # Compute eigenvalues of the sample covariance matrix
  ev <- svd(Y)$d^2 / N  
  N <- length(ev)

  if (is.null(rmax))                    # Set a value for rmax
    rmax <- 3 * sqrt(N)
  j <- rmax + 1

  diffs <- ev - c(ev[-1], 0)            # adjacent eigenvalues difference

  for (i in 1:niter) {
    y <- ev[j:(j+4)]
    x <- ((j-1):(j+3))^(2/3)
    lm.coef <- lm(y ~ x)
    delta <- 2 * abs(lm.coef$coef[2])   # choose delta in step 2

    # select the the adjacent eigen differences larger than delta
    idx <- which(diffs[1:rmax] > delta)

    if (length(idx) == 0)
      hatr <- 0
    else hatr <- max(idx) # select the largest eigen difference

    newj = hatr + 1

    if (newj == j) break
    j = newj              # Repeat until convergence  
  }
  return(hatr)
}
{{< /highlight >}}

## 3. Strong factor estimation and robustness check


Simulation design: we use the <a href="#DGP_strong"> same DGP </a> for generating strong factors as the IC method:

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

The corresponding R code of applying the ED method for this DGP is:

{{<highlight r>}}

S <- 10                   # Number of Simulation
N_set <- c(50, 100)
T_set <- c(50, 100, 200)
NT_comb <- expand.grid(N_set, T_set)  # Combination of N and T pairs

r_hat_ED <- matrix(NA, S, 4)
Factor_ED <- matrix(NA,nrow(NT_comb),6)
colnames(Factor_ED) <- c("N","T", "white", "serial", "cross", "gamma")

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

    # Apply the ED method to estimate the number of strong factors
    # in the DGP:
    r_hat_ED[s,1] <-  EigenDiff(X, niter=10, rmax= 15)
    r_hat_ED[s,2] <-  EigenDiff(X1, niter=10, rmax= 15)
    r_hat_ED[s,3] <-  EigenDiff(X2, niter=10, rmax= 15)
    r_hat_ED[s,4] <-  EigenDiff(X3, niter=10, rmax= 15)
  }
    Factor_ED[i,] <- c(N,T,colMeans(r_hat_ED))
}
{{< /highlight >}}

The results of the ER method for estimating the number of strong factors with different types of error terms are:

{{<highlight r>}}
       N   T white serial cross gamma
[1,]  50  50   5.1    5.0   4.3   4.7
[2,] 100  50   5.0    5.4   5.1   6.0
[3,]  50 100   5.0    4.9   4.9   5.9
[4,] 100 100   5.0    5.0   5.0   6.1
[5,]  50 200   5.0    5.0   5.0   6.1
[6,] 100 200   5.0    5.0   5.0   6.1
{{< /highlight >}}

From the results above, we can see that the ED method is robust when the error terms in the factor model are high serially and cross-sectionally correlated. It is not robust when we have non-gaussian distributed error terms in the factor model, however.

## 4. Weak factor estimation

An advantage of the ED method comparing with the IC and ER methods is that the consistency of the ED estimator can allow for much weaker strength of the factors. To show this, let's apply the ED method to <a href="#DGP_weak"> the DGP </a> that we used to generate weak factors before:

{{<highlight r>}}
Factor_ED_weak <- matrix(NA,nrow(NT_comb),3)
colnames(Factor_ED_weak) <- c("N","T","r_hat_NE")

for(i in 1:nrow(NT_comb)){
  for(s in 1:S){

    # Data generating process
    N <- NT_comb[i,1]
    T <- NT_comb[i,2]
    X <- DGP_weak(N,T,5)$Y

    # Apply the ED method to estimate the number of
    # strong factors in the DGP:    
    r_hat_ED[s,] <-  EigenDiff(X, niter=10, rmax= 15)
  }
    Factor_ED_weak[i,] <- c(N,T,mean(r_hat_ED))
}
{{< /highlight >}}

The goal is to estimation the number of strong and useful weak factors in the DGP. The results of the ED method to estimate the number of weak factors in the DGP are:

{{<highlight r>}}
       N   T r_hat_ED
[1,]  50  50      2.7
[2,] 100  50      3.5
[3,]  50 100      3.0
[4,] 100 100      3.4
[5,]  50 200      3.3
[6,] 100 200      3.5
{{< /highlight >}}

From the results above, we can see that the ED method can estimate the number of strong and useful weak factors in the DGP very well in finte samples.
