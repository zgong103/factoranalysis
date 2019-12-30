---
title: "Strong and weak factor assumptions"
date: 2019-12-27
tags: ["assumptions", "introduction"]
draft: false
---

## Strong and weak factor assumptions

In this section, we briefly go over the strong and weak factor assumptions. Assuming factors $F_{t}$ and noise $e_{t}$ are uncorrelated and have zero mean, and normalization $\mathbb{E}(F_{t}F_{t}^{'}) = I_{r}$ for identification, then the population covariance matrix of the factor model (1) can be expressed as
\begin{equation}
    \Sigma_{Y} = LL^{'} + \Sigma_{e}, \quad \quad \quad (2)
\end{equation}
where $\Sigma_{Y}$ and $\Sigma_{e}$ are the $N \times N$ population covariance matrix of $Y_{t}$ and $e_{t}$, respectively.

**Strong Factor assumption:**

`For (2), we assumed that $L^{'}L/N \rightarrow \Sigma_{L}$ for some $r_{0} \times r_{0}$ positive definite matrices $\Sigma_{L}$ and all the eigenvalues of $\Sigma_{e}$ are bounded as $N,T \rightarrow \infty$.`

Under this assumption, the top $r_{0}$ eigenvalues of $\Sigma_{Y}$ are diverge at the rate $O(N)$ while the rest of its eigenvalues are bounded as $N, T \rightarrow \infty$. This is the critical assumption for those methods to consistently estimate the number of strong factors as $N, T \rightarrow \infty$.

**Weak Factor assumption:**

`In contrast to strong factors, for the weak factors, we assumed that $L^{'}L \rightarrow \Sigma_{L}$ instead of $L^{'}L/N \rightarrow \Sigma_{L}$ and all the eigenvalues of $\Sigma_{e}$ are bounded as $N,T \rightarrow \infty$.`

Under this assumption, all the eigenvalues of $\Sigma_{Y}$ are bounded as $N,T \rightarrow \infty$ and PCA or MLE estimators for estimating factors and corresponding loadings are not consistent.
