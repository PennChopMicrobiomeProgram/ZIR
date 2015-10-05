#' Modified Wilcoxon rank test (ZIW) for zero-inflated data
#'
#' @param x a vector of data for one group
#' @param y a vector of data for another group
#' @param perm use permutation to calculate the pvalue
#' @param perm.n number of permutations used for pvalue calculation
#' @return modified wilcoxon rank sum statistic and pvalue
#' @export
#' @import stats utils graphics methods
#' @examples
#' x <- c(rep(0,5),rlnorm(20, meanlog = 0, sdlog = 1))
#' y <- c(rep(0,10),rlnorm(20, meanlog = 2, sdlog = 1))
#' ziw(x, y, perm = FALSE)
#' ## use permutations to calculate the pvalue
#' ziw(x, y, perm = TRUE)
#' @references Wanjie Wang, Eric Z. Chen and Hongzhe Li (2015). Rank-based tests for compositional distributions with a clump of zeros. Submitted.

ziw = function(x, y, perm = FALSE, perm.n = 10000) {
  calculate_ziw_statistic = function(x,y){
    ## total observations in each vector
    N <- c(length(x), length(y))
    ## number of non-zero observations in each vector
    n <- c(sum(x != 0), sum(y != 0))
    ## non-zero proportion
    prop <- n / N 
    pmax <- max(prop) 
    pmean <- mean(prop)
    ## keep only round(pmax * N) observations in each group
    Ntrun <- round(pmax * N)
    ## truncate zeros in X and Y
    Xtrun <- c(x[x != 0], rep(0, Ntrun[1] - n[1]))
    Ytrun <- c(y[y != 0], rep(0, Ntrun[2] - n[2]))
    rankdata <- sum(Ntrun) + 1 - rank(c(Xtrun, Ytrun))
    ## mean of the modified wilcoxon rank sum statistic
    r <- sum(rankdata[1:Ntrun[1]])
    s <- r - Ntrun[1] * (sum(Ntrun) + 1) / 2
    ## variance of the modified wilcoxon rank sum statistic
    v1   <- pmean * (1 - pmean)
    var1 <- N[1] * N[2] * pmean * v1 * (sum(N) * pmean + 3 * (1 - pmean) / 2) + 
      sum(N^2) * v1^2 * 5 / 4 + 2 * sqrt(2 / pi) * pmean * (sum(N) * v1) ^ (1.5) *
      sqrt(N[1] * N[2])
    var1 <- var1 / 4
    var2 <- N[1]*N[2]*pmean^2*(sum(N)*pmean+1)/12
    vars <- var1 + var2
    ## modified wilcoxon rank sum statistic
    w <- s / sqrt(vars)
    return(w)
  }
  
  ## calclulate the ziw statistic
  w <- calculate_ziw_statistic(x,y)
  
  ## calculate the pvalue
  if (perm) {
    numrep <- 10000 
    permu.w <- rep(0, numrep)
    N <- c(length(x), length(y))
    Z <- c(x, y)
    for (i in 1:numrep) {
      set.seed(i)
      ind <- sample(sum(N), N[1])
      x <- Z[ind] 
      y <- Z[-ind]
      permu.w[i] <- calculate_ziw_statistic(x,y) 
    }
    p <- sum(abs(w) < abs(permu.w)) / numrep
  }
  else{
    p <- 2 * (1 - pnorm(abs(w)))
  }
  
  return(list(p.value = p, statistics = w))
}
