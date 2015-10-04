#' Modified Kruskal Wallis test (ZIKW) for zero-inflated data
#'
#' @param x a numeric vector of data values with groups defined by the group parameter, or a list of numeric data vectors with each vector in the list as one group.
#' @param group a vector giving the groups for the corresponding elements of x. Ignored if x is a list.
#' @param perm use permutations to calculate the pvalue
#' @param perm.n number of permutations used for pvalue calculation
#' @return modified Kruskal Wallis test statistic and pvalue
#' @export
#' @examples
#' ## x is a list
#' x <- list(group1 = c(rep(0,5),rlnorm(20, meanlog = 0, sdlog = 1)),
#'      group2=c(rep(0,10),rlnorm(20, meanlog = 1, sdlog = 1)),
#'      group3=c(rep(0,15),rlnorm(20, meanlog = 2, sdlog = 1)))
#' zikw(x, perm = FALSE)
#' ## x is a vector
#' x <- c(c(rep(0,5),rlnorm(20, meanlog = 0, sdlog = 1)),
#'     c(rep(0,10),rlnorm(20, meanlog = 1, sdlog = 1)),
#'     c(rep(0,15),rlnorm(20, meanlog = 2, sdlog = 1)))
#' group <- c(rep('group1',25),rep('group2',30),rep('group3',35))
#' zikw(x, group, perm = FALSE)
#' ## use permutations to calculate the pvalue
#' \dontrun{zikw(x, group, perm = TRUE)}

zikw <- function(x, group, perm = FALSE, perm.n = 10000) {
  ## transform x into a list
  vect2list = function(x,group){
    s <- unique(group)
    x.list <- list()
    for (i in s) {
      x.list[[i]] <- x[group == i]
    }
    return(x.list)
  }
  ## check if x is a list
  if (class(x) == "numeric" || class(x) == "data.frame") {
    x <- vect2list(x,group)
  }
  
  ## this function takes a list as the input
  calculate_zikw_statistic = function(x){
    ## number of groups
    K <- length(x)
    ## total observations in each group
    N <- rep(0, K)
    ## number of non-zero observations in each group
    n <- rep(0, K)
    xvec <- numeric(0)
    ## count total and non-zero observations in each group
    for (i in 1:K) {
      N[i] <- length(x[[i]])
      n[i] <- sum(x[[i]] != 0)
      xvec <- c(xvec, x[[i]])
    }
    ## non-zero proportion
    prop <- n / N
    pmax <- max(prop)
    ## keep only round(pmax * N) observations in each group
    Ntrun <- round(pmax * N)
    ## truncate zeros in each group
    Xtrun.vec <- numeric(0)
    for (i in 1:K) {
      data <- x[[i]]
      Xtrun.vec <-
        c(Xtrun.vec, data[data != 0], rep(0, Ntrun[i] - n[i]))
    }
    rankdata <- sum(Ntrun) + 1 - rank(Xtrun.vec)
    r <- sum(rankdata[1:Ntrun[i]])
    
    for (i in 2:K) {
      r <- c(r, sum(rankdata[1:Ntrun[i] + sum(Ntrun[1:(i - 1)])]))
    }
    s <- r - Ntrun * (sum(Ntrun) + 1) / 2
    u <- numeric(0)
    
    for (i in 1:(K - 1))
      u <- c(u, N[i + 1] * sum(s[1:i]) - sum(N[1:i]) * s[i + 1])
    u <- u / sum(N) ^ 2
    thetam <- mean(prop)
    simun <- matrix(0, nrow = 5000, ncol = K)
    simup <- simun
    for (ss in 1:K) {
      simun[,ss] <- rbinom(5000, N[ss], thetam)
      simup[,ss] <- simun[,ss] / N[ss]
    }
    simupmax <- apply(simup, 1, max)
    varsimu <- numeric(K - 1)
    
    varsimu[1] <-
      N[2] ^ 2 * mean(simupmax ^ 2 * (simup[,1] - simup[,2]) ^ 2) * N[1] ^ 2
    varu2 <- N[2] * N[1] * (N[1] + N[2])
    
    for (ss in 2:(K - 1)) {
      varsimu[ss] <-
        N[ss + 1] ^ 2 * mean(simupmax ^ 2 * (apply(simun[,1:ss], 1, sum) - simup[,ss +
                                                                                   1] * sum(N[1:ss])) ^ 2)
      varu2 <-
        c(varu2, N[ss + 1] * sum(N[1:ss]) * sum(N[1:(ss + 1)]))
    }
    varsimu <- varsimu / (sum(N)) ^ 2 / 4
    
    varu2 <-
      varu2 * thetam ^ 2 * (thetam + 1 / sum(N)) / 12 / (sum(N)) ^ 2
    varu <- varsimu + varu2
    ## modified Kruskal Wallis test statistic
    w <- sum(u ^ 2 / varu)
    return(w)
  }
  
  ## calclulate the zikw statistic
  w <- calculate_zikw_statistic(x)
  
  ## calculate the pvalue
  if (perm) {
    numrep <- perm.n
    permu.w <- rep(0, numrep)
    perm.group <- rep(c(1:length(x)),sapply(x,length))
    xvec <- unlist(x)
    ## need to permute a list
    for (i in 1:numrep) {
      set.seed(i)
      xvec.perm <- sample(xvec, length(xvec))
      x.list.perm <- vect2list(xvec.perm,group)
      permu.w[i] <- calculate_zikw_statistic(x.list.perm)
    }
    pw <- sum(abs(w) < abs(permu.w)) / numrep
  }
  else{
    K <- length(x)
    pw <- pchisq(w, K - 1, lower.tail = F)
  }
  
  return(list(p.value = pw, statistics = w))
}
