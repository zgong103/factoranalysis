---
title: "ER Method"
date: 2019-12-29
tags: ["strong", "methods"]
draft: false
---

## 1. Introduction

The ER estimator is defined as

 `$$ \hat{r}_{\mathrm{ER}}=\operatorname{argmin}_{0 \leq r \leq r_{max}} \lambda_{r}/\lambda_{r+1},$$`
with $\lambda_{0} = \sum_{r=1}^{\min(N,T)}\lambda_{r}/ \log \min(N,T)$.

The intuition for this method to work is very simple: based on strong factor assumption, for any $j \neq r_{0}$ the ratio $\lambda_{j}/ \lambda_{j+1}$ converges to $O(1)$ as $N, T \rightarrow \infty$, while the the ratio $\lambda_{r_{0}}/\lambda_{r_{0}+1}$  diverges to infinity. For the details of this method, please refer [Ahn and Horenstein, 2013](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA8968).

## 2. R codes

Based on those forms, we can write the R code for the ER method as follow:

{{< highlight r >}}

EigenRatio <- function(X){

  N <- length(X[,1]);
  T <- length(X[1,]);

  if(N<=n){
    S = X%*%t(X)/T;
  } else {
    S = t(X)%*%X/N;
  }

  # Generate eigenvalues for sample covariance matrix
  eig <- sort(eigen(S)$value, decreasing = TRUE)

  m <- min(T, N)

  # lamda=0 case
  lamda_0 <- sum(eig)/log(m)
  eig_plus <- c(lamda_0, eig)

  # calculate the value of r_max
  dif <- eig - (sum(eig[1:m])/m)
  plus <- sum(as.numeric(dif >= 0))
  r_max_ER <- min(plus, 0.1*m)

  # generate eigen ratios for the sequence of r values
  eigen_ratio <- eig_plus[1:(r_max_ER+1)]/eig_plus[2:(r_max_ER+2)]

  # select the largest eigen ratios
  k_hat <- which.max(eigen_ratio) - 1
  return(k_hat)
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

The corresponding R code of applying the ER method for this DGP is:

{{<highlight r>}}

S <- 10                   # Number of Simulation
N_set <- c(50, 100)
T_set <- c(50, 100, 200)
NT_comb <- expand.grid(N_set, T_set)  # Combination of N and T pairs

r_hat_ER <- matrix(NA, S, 4)
Factor_ER <- matrix(NA,nrow(NT_comb),6)
colnames(Factor_ER) <- c("N","T", "white", "serial", "cross", "gamma")

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

    # Apply the ED method to estimate the number of
    # strong factors in the DGP:

    r_hat_ER[s,1] <-  EigenRatio(X)
    r_hat_ER[s,2] <-  EigenRatio(X1)
    r_hat_ER[s,3] <-  EigenRatio(X2)
    r_hat_ER[s,4] <-  EigenRatio(X3)
  }
    Factor_ER[i,] <- c(N,T,colMeans(r_hat_ER))
}
{{< /highlight >}}

The results of the ER method for estimating the number of strong factors with different types of error terms are:

{{<highlight r>}}
       N   T white serial cross gamma
[1,]  50  50     5    5.0     5   1.8
[2,] 100  50     5    4.9     5   2.2
[3,]  50 100     5    5.0     5   1.3
[4,] 100 100     5    5.0     5   6.0
[5,]  50 200     5    5.0     5   0.7
[6,] 100 200     5    5.0     5   6.0
{{< /highlight >}}

From the results above, we can see that the ER method is quite robust when the error terms in the factor model are high serially and cross-sectionally correlated. It is not robust when we have non-gaussian distributed error terms in the factor model, however.

## 4. Weak factor estimation

Besides, the ER method is designed to estimation the number of factors under strong factor assumption, so it may fail to detect weak factors in the data. To show this, let's apply the IC method to <a href="#DGP_weak"> the DGP that we used to generate weak factors </a>:


{{<highlight r>}}
Factor_ER_weak <- matrix(NA,nrow(NT_comb),3)
colnames(Factor_ER_weak) <- c("n","T","r_hat_ER")

for(i in 1:nrow(NT_comb)){
  for(s in 1:S){

    # Data generating process
    N <- NT_comb[i,1]
    T <- NT_comb[i,2]
    X <- DGP_weak(N,T,5)$Y

    # Apply the ED method to estimate the number of strong
    # factors in the DGP:
    r_hat_ER[s,] <-  EigenRatio(X)
  }
    Factor_ER_weak[i,] <- c(N,T,mean(r_hat_ER))
}
{{< /highlight >}}

The goal is to estimation the number of strong and useful weak factors in the DGP. The results of the ER method to estimate the number of weak factors in the DGP are:

{{<highlight r>}}
       N   T r_hat_ER
[1,]  50  50        1
[2,] 100  50        1
[3,]  50 100        1
[4,] 100 100        1
[5,]  50 200        1
[6,] 100 200        1
{{< /highlight >}}

From the results above, we can see that the ER method fails to detect the number of weak factors in the DGP in finite samples, but it can precisely seprate the number of strong factors from weak ones in the DGP.
