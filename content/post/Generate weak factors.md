---
title: "Generate weak factors"
date: 2019-12-28
tags: ["weak", "DGP"]
draft: false
---

Before introducing the methods to estimate the number of strong and weak factors for high dimensional data, let me first introduce how to <a id="DGP_weak"> genernate the weak factors </a> in the DGP and corresponding R codes. Consider the factor model as:

$$
Y = X + \Sigma^{\frac{1}{2}}E = \Sigma^{\frac{1}{2}} (\sqrt{T} \hat{U} \hat{D} \hat{V}' + E),
$$

where $\sqrt{T} \hat{U} \hat{D} \hat{V}^{'}$ is the singular value decomposition (SVD) for `$\Sigma^{-\frac{1}{2}}X$ with $\hat{U} \in \mathbb{R}^{N\times \min (N,T)}$, $\hat{V} \in \mathbb{R}^{T \times \min (N,T)}$,  $\hat{D}= \text{diag} (\hat{d}_{1}, \hat{d}_{2}, \cdots, \hat{d}_{\min (N, T)})$, $\hat{U}' \hat{U} = \hat{V}' \hat{V} = I_{\min (N, T)}$, and $\hat{d}_{1} \geq \hat{d}_{2} \geq , \cdots, \geq \hat{d}_{\min (N, T)}$.`

Using the guidance from RMT, we can generate the factors with different strength by specified the entries in $\hat{D}$. Specifically, in our DGP as an example, we generate one strong factor, two useful weak factors, one harmful weak factor and one undetectable weak factor by following the thresholds below:

* **Undetectable factor**, `$d^{2}_{i} < \mu_{F}$`, the factor is asymptotically undetectable by PCA or MLE based methods.

* **Harmful weak factor**, `$\mu_{F} < d^{2}_{i} < \mu_{F}^{*}$`, including the factor in the model will make the loss for recovering signal matrix $X$ larger.

* **Useful weak factor**, `$\mu_{F}^{*} < d^{2}_{i} = O(1)$`, including the factor will reduce the loss for recovering signal matrix $X$.

* **Strong factor**, `$d^{2}_{i}$` grows proportionally to $N$.

where $\mu_{F}=\sqrt{\gamma}$ and
$$\mu_{F}^{\star}=\frac{1+\gamma}{2}+\sqrt{\left(\frac{1+\gamma}{2}\right)^{2}+3 \gamma}, \quad  \text{for} \quad \frac{N}{T} \rightarrow \gamma$$

are detection and estimation thresholds based on the results from RMT (Random matrix theory). For the error term $\Sigma^{\frac{1}{2}}E$ in the factor model, we assume homoscedastical noise `$\Sigma= I_{N}$` and `$E=\left(e_{i t}\right)_{N \times T}: e_{it} \stackrel{\text { iid }}{\sim} \mathcal{N}(0,1)$` for simplicity. Based on the factor model and the thresholds to generate weak factors, we can write the R code for our DGP as follow:

{{<highlight r>}}
# Generate data
# Imputs: N is the dimension, n is the sample size,
# k is the number of factors
rm(list = ls())
set.seed(111)

DGP_weak <-  function (N, T, k) {

  r <- T/N;


  u_star <- (1/(sqrt(r)));      # compute the dectection threshold         

  u_starF <- (1+1/r)/2+sqrt(((1+1/r)/2)^2 + 3/r);      # compute the estimation threshold

#  print("Detection threshod, estimating threshold:")
#  print( c(u_star, u_starF))

  Sigma <- rep(1, N)                              # white noise case
  noise <- matrix(rnorm(N * T), nrow = N)         # Generate the noise

  # Generate orthogonal matrix to construct signal matrix
  U <- randOrtho(N,k)                      
  V <- randOrtho(T,k)

  # Generate factors with different strength
  D <- c(u_star/2,(-u_star+u_starF)/2 + u_star,
         u_starF*c(1.5, 2.5), N * 1.5)
  D <- diag(sort(D, decreasing = TRUE))

  # Construct the signal matrix X
  new.X <- sqrt(Sigma)^{-1} * U %*% sqrt(D) %*% t(V);
  U1 <- svd(new.X, nu = k, nv = k)$u;
  adjust.X <- sqrt(Sigma) * U1 %*% sqrt(D) %*% t(V)

  Y <- sqrt(T) * adjust.X + noise;                # DGP for weak factors

  # SVD for reweighted signal matrix
  svd.X <- svd(adjust.X, nu = k, nv = k);         
  U <- svd.X$u;
  V <- svd.X$v;
  D <- diag(svd.X$d[1:k]^2);

  return (list(Y=Y, D = sqrt(D), U=U, V=V))
}
{{< /highlight >}}

Function for generating the orthogonal matrix $U$ and $V$:

{{<highlight r>}}
randOrtho <- function(N, k) {
  Z <- matrix(rnorm(N * k), nrow=N);
  Q1 <- qr.Q(qr(Z))[, 1:k];            # QR decomposition
  if (k==1) {
    S <- sample(c(1, -1), k, replace=T)
    return(S*Q1)
  }
  else {
    S <- diag(sample(c(1,-1),k,replace=T));
    return(Q1%*%S)
  }
}
{{< /highlight >}}

Example: we can use the R codes above to generate a signal matrix $Y$ with $N=T=50$ and five factors, in which we have one strong factor, two useful weak factors, one harmful weak factor and one undetectable weak factor:

{{<highlight r>}}
Y <- DGP_weak(50, 50, 5)

(Y$D)^2

     [,1] [,2] [,3] [,4] [,5]
[1,]   75  0.0  0.0    0  0.0
[2,]    0  7.5  0.0    0  0.0
[3,]    0  0.0  4.5    0  0.0
[4,]    0  0.0  0.0    2  0.0
[5,]    0  0.0  0.0    0  0.5
{{< /highlight >}}
