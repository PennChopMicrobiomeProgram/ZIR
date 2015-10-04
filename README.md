# ZIR (Zero-Inflated Rank Test)

## Introduction
When the underlying distribution of data is unknown, usually we use nonparametric test to test the difference between groups. For example, Wilcoxon rank sum test (Mann-Whitney U test) is often used for two group comparison, and Kruskal Wallis test is used for multiple group comparison. However, sometimes the data may contain many zeroes and the traditional rank based test, such as Wilcoxon rank sum test and Kruskal Wallis test have low power due to the ties in the data. 

Here we developed a R package called ZIR (Zero-Inflated Rank Test), which has more power than the traditional rank based test when the data contain many zeroes, and have similar power as the traditional rank based test when the data don't have too many zeros.

## Installation
You can install our ZIR package from Github
```r
install.packages("devtools")
devtools::install_github("chvlyl/ZIR")
library(ZIR)
```

## Basic Usage
### Modified Wilcoxon rank sum test (ZIW) for zero-inflated data
```r
x <- c(rep(0,5),rlnorm(20, meanlog = 0, sdlog = 1))
y <- c(rep(0,10),rlnorm(20, meanlog = 2, sdlog = 1))
ziw(x, y, perm = FALSE)
```
If you want to use permutations to generate pvalues, you can set perm to TRUE
```r
ziw(x, y, perm = TRUE)
```

### Modified Kruskal Wallis test (ZIKW) for zero-inflated data
```r
x <- list(group1 = c(rep(0,5),rlnorm(20, meanlog = 0, sdlog = 1)),
      group2=c(rep(0,10),rlnorm(20, meanlog = 1, sdlog = 1)),
      group3=c(rep(0,15),rlnorm(20, meanlog = 2, sdlog = 1)))
zikw(x, perm = FALSE)
```
If you want to use permutations to generate pvalues, you can set perm to TRUE
```r
zikw(x, perm = TRUE)
```

## Citation
Wanjie Wang, Eric Z. Chen and Hongzhe Li (2015). Rank-based tests for compositional distributions with a clump of zeros. Submitted.