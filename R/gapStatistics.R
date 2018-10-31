##' gap statistics to tune gamma
##'
##' gap statistics to tune gamma, the tuning parameter to control number of features.
##' @title gap statistics
##' @param d A list of S studies, each study is a combined data matrix n*J, where n is number of subjects, J=J1+J2+... and J1 is number of features in omics dataset 1 and J2 is number of features in omics dataset 2...
##' @param K number of clusters
##' @param B number of permutations.
##' @param gamma Penalty on total number of features. Default is given. Larger gamma will yeild small number of selected features.
##' @param alpha balance between group sparsity and individual sparsity. alpha=1 yeilds no group sparsity. alpha=0 yeilds no individual penalty.
##' @param group Prior group information. Potentially these group can contain overlap features. group is a list and each element of the list is feature index.
##' @param seed random seed
##' @param silence do not print progress in silence = TRUE.
##' @return a table. Each row represents a gamma. 
##' \item{gamma}{input gammas}
##' \item{score}{objective score: see obj0 in ISKmeans function}
##' \item{E.score}{mean value of permutated score}
##' \item{se.score}{standard error of permutated score}
##' \item{gapStat}{gap statistics, score - E.score}
##' \item{numF}{number of selected features}
##' @export
##' @author Caleb
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
##' gapResult <- gapStatistics(d=S,K=3,B=10,group=group)
##' print(gapResult)
gapStatistics <-
function(d,K=3,B=10,gamma=NULL,alpha=1, group=NULL,seed=15213,silence=FALSE){
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0)
     stop("'B' has to be a positive integer")
  if (is.null(gamma))
	  gamma <- seq(0,0.8, 0.05)
  if (min(gamma) < 0)
      stop("gamma should be greater than or equal to 0")

  ## get true objective score
  if (!silence)
  	cat("calculating true score...\n")
  set.seed(seed)
  trueRes <- MISKmeans(d,K=K,gamma=gamma,alpha=alpha,group=group,silent=TRUE)

  numF <- sapply(trueRes,function(x) sum(x$ws!=0))
  nGamma <- sapply(trueRes,function(x) x$groupInfo$gamma)
  score <- sapply(trueRes,function(x) x$obj0)
  groupInfos <- lapply(trueRes,function(x) x$groupInfo)
  
  if (!silence)
     cat("calculating permutated score, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", sep = "")
   E.score.full <- NULL
  for(b in 1:B){
	  cat(".", if (b%%50 == 0)
	         paste(b, "\n"))

	  set.seed(15213+b)
      ad = lapply(d,function(dd) t(apply(dd,1,function(x) sample(x))))

	  permRes <- MISKmeans(ad,K=K,gamma=gamma,alpha=alpha,group=group,penaltyInfo=groupInfos,silent=TRUE)
	  scoreperm <- sapply(permRes,function(x) x$obj0)
	  E.score.full <- rbind(E.score.full, scoreperm)
  }

  if (!silence && (B%%50 != 0))
          cat("", B, "\n")

  E.score <- colMeans(E.score.full)
  se.score <- apply(E.score.full,2,sd)/nrow(E.score.full)

  gapStat <- score - E.score

  stat <- data.frame(gamma, score, E.score, se.score, gapStat,numF=numF)
  bestGamma <- gamma[which.min(gapStat)]
  result <- list(stat=stat, bestGamma=bestGamma)
  return(result)
}
