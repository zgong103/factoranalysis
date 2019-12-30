---
title: "Basic Factor Model"
date: 2019-12-26
tags: [ "introduction", "factor model"]
draft: false
---

## Basic factor model

Factor analysis based on a model that separates the observed data into an unobserved systematic part (signal part) and an unobserved error part (noise part). The systematic part captures the main information of the data so that we want to separate it from noise part. Specifically, let $Y_{it}$ be the observed data for the $i$-th cross-section unit at time $t$, for $i=1,2,\cdots, N$ and $t=1,\cdots,T$. The factor model for $Y_{it}$ is given by
\begin{equation}
    Y_{it} =  L_{i}^{'} F_{t} + e_{it},  
\; i = 1,\dots,N, t = 1,\dots, T, \quad \quad \quad (1)
\end{equation}
 where $F_{t}$ is a $(r_{0} \times 1)$ vector of common factors, $L_{i}$ is a $(r_{0} \times 1)$ vector of loadings associated with $F_{t}$, and $e_{it}$ is the idiosyncratic component (noise part) of $Y_{it}$. The number of true factors in the model is $r_{0}$. The product of $L_{i}^{'}F_{t}$ is called the common component (signal part) of $Y_{it}$. The factors, their loadings, as well as the idiosyncratic errors are not observable. The goal is to estimate the number of factors $r_{0}$ in the observed data under strong and weak factor assumptions.
