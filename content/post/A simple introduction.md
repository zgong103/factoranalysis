---
title: "A Simple Introduction"
date: 2019-12-24
tags: ["introduction"]
draft: false
---

## Introduction

In the following posts, we provide a systematic review and explanations for the R codes that used in current popular methods for choosing the number of strong and weak factors in high dimensional data. Examples for how to apply those methods in practice are included. Specifically, those methods we reviewed are lists as follow:

* Methods for estmating the number of strong factors:
   + [IC Method: information criteria based methods (Bai and Ng, 2002).](http://zgong103.github.io/factor.github.io/post/goisforlovers/)
   + [ER Method: eigenvalue ratio based method (Ahn and Horenstein, 2013).](#ER)
   + [ED Method: eigenvalue difference based method (Onatski, 2010).](#ED)

* Methods for estmating the number of weak factors:
  + [NE Method: sample eigenvalue based method (Nadakuditi and Edelman, 2008).](#NE)
  + [BCV Method: bi-cross-section based method (Owen and Wang, 2015).](#BCV)

Most methods for choosing the number of factors are based on the results from random matrix theory(RMT), which studies the distribution of sample eigenvalues and requires i.i.d and gaussian assumption on the error terms. These restrictions may not appropriate when we want to apply them in practice. In this page, we also show that those methods we reviewed are not robust in the simulation when the error terms in the factor model are serially and cross-sectionally correlated or have non-gaussian distributions.
