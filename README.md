# ZIR (Zero-Inflated Rank Test)

## Desciprtion
When the underlying distribution of data is unknown, usually we use nonparametric test to test the differece between groups. For example, Wilcoxon rank sum test (Mann-Whitney U test) is often used for two group comparison, and Kruskal Wallis test is used for multiple group comparison. However, sometimes the data may contain many zeroes and the tranditional rank based test, such as Wilcoxon rank sum test and Kruskal Wallis test have low power due to the ties in the data. 

Here we developed a R package called ZIR (Zero-Inflated Rank Test), which has more power than the tranditional rank based test when the data contain many zeroes, and have similar power as the tranditional rank based test when the data don't have too many zeros.

## Install
You install them from CRAN with install.packages("x")

You can install our ZIR package from Github
```r
install.packages("devtools")
devtools::install_github("chvlyl/ZIR")
```

## Use