##' MISKmenas
##'
##' Estimate a propriate tuning parameter alpha
##' @title Infer alpha
##' @param resMISKmeans The result returned from MISKmeans.
##' @param balance selecting alpha such that the l1 norm penalty = balance * group penalty
##' @return alpha values
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
##' inferAlpha(res) 
##' 	
inferAlpha <- function(resMISKmeans, balance=1){
	groupInfo <- resMISKmeans$groupInfo
	ws <- resMISKmeans$ws
	alpha <- groupInfo$alpha
	groupLevel <- groupInfo$groupLevel
	split_pos <- split(groupInfo$genePos, groupLevel)
	split_coef <- split(groupInfo$coef, groupLevel)
	allPenalty_perGroup <- Map(function(x,y){
		sqrt(sum((ws[x] * y)^2))
	}, split_pos, split_coef			
		)
	allPenalty_combined <- Reduce("+", allPenalty_perGroup)
	l1Penalty <- groupInfo$alpha * groupInfo$gamma * sum(abs(ws))
	groupPenalty <- allPenalty_combined - l1Penalty
	l1Penalty_2 <- l1Penalty/alpha
	groupPenalty_2 <- groupPenalty/(1-alpha)
			
	alphas <- groupPenalty_2*balance/(groupPenalty_2*balance + l1Penalty_2)
	alphas
}
	