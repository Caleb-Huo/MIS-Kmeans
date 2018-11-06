##' MISKmenas
##'
##' Estimate a propriate tuning parameter lambda
##' @title Infer lambda
##' @param resMISKmeans The result returned from MISKmeans.
##' @param balance selecting lambda such that the separation ability = balance * matching ability
##' @return lambda values
##' @author Caleb
##' @export
##' @examples
##' 
##' S <- 2
##' K <- 3
##' G <- 1000
##' g1 <- 50
##' g2 <- 50
##' n0 <- 20
##' n <- K*n0
##' labels <- cut(1:n,breaks=K,labels=FALSE)
##' 
##' set.seed(32611)
##' S1 <- matrix(rnorm(G*n), nrow=G, ncol=n)
##' S2 <- matrix(rnorm(G*n), nrow=G, ncol=n)
##' 
##' S1[1:g1, labels==1] <- S1[1:g1, labels==1] + 2
##' S1[1:g1, labels==3] <- S1[1:g1, labels==3] - 2
##' S1[g1 + 1:g2, labels==1] <- S1[g1 + 1:g2, labels==1] - 2
##' S1[g1 + 1:g2, labels==2] <- S1[g1 + 1:g2, labels==2] + 2
##' 
##' S2[1:g1, labels==2] <- S2[1:g1, labels==2] + 2
##' S2[1:g1, labels==1] <- S2[1:g1, labels==1] - 2
##' S2[g1 + 1:g2, labels==2] <- S2[g1 + 1:g2, labels==2] - 2
##' S2[g1 + 1:g2, labels==3] <- S2[g1 + 1:g2, labels==3] + 2
##' 
##' S = list(t(S1),t(S2))
##' groups <- Map('c',1:g1,g1 + 1:g2)
##' 
##' res <- MISKmeans(d = S, K = 3, gamma = 0.4, group = groups)
##' 
##' inferLambda(res) 
##' 
inferLambda <- function(resMISKmeans, balance=1/2){
	lambdas <- sum(resMISKmeans$ws * resMISKmeans$per_ratio)/sum(resMISKmeans$ws * resMISKmeans$per_match)*balance
	lambdas
}
