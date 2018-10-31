##' MISKmenas
##'
##' Meta-analytic multi-level omics data integration with incorporation of prior group structure
##' @title MISKmeans
##' @param d combined data matrix n*J, where n is number of subjects, J=J1+J2+... and J1 is number of features in omics dataset 1 and J2 is number of features in omics dataset 2...
##' @param K number of clusters
##' @param gamma Penalty on total number of features. Larger gamma will yeild small number of selected features.
##' @param lambda A tuning parameter controlling the balance between separation ability (BCSS/TSS) and matching function. lambda is set to be 0.5 by default.
##' @param alpha balance between group sparsity and individual sparsity. alpha=1 yeilds no group sparsity. alpha=0 yeilds no individual penalty. Default alpha=0.5
##' @param group Prior group information. Potentially these group can contain overlap features. group is a list and each element of the list is feature index.
##' @param nstart Number of initials for Kmeans for sparse Kmeans
##' @param wsPre Initial feature weight.
##' @param penaltyInfo only for the purpose of gap statitics. Here we will fix the penalty design to perform gap statistics. The input should be a list of groupInfo. See groupInfo for details.
##' @param silent Output progress.
##' @param maxiter Maximum numbre of iteration between ws and Cs.
##' @param sampleSizeAdjust logical argument,
##' controlling whether to adjust for sample size.
##' If true, that means study with larger sample size will have a larger impact.
##' If false, each study has equal contribution.
##' Without prior information, sampleSizeAdjust=FALSE is suggested since we are not sure about data quality.
##' @return m lists, m is length of gamma parameters. If gamma is a scalar, m = 1. The following items are included in the list.
##' \item{ws}{weight for each feature. Zero weight means the feature is not selected.}
##' \item{Cs}{Cluster Assignment}
##' \item{obj0}{sum of weighted separation ability. This term is for the purpose of gap statistics.}
##' \item{objective}{objective value}
##' \item{groupInfo}{a list containing group design, alpha, gamma}
##' @author Caleb
##' @useDynLib MISKmeans
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
##' heatmap(S1, Rowv=NA, Colv=NA)
##' 
##' S2[1:g1, labels==2] <- S2[1:g1, labels==2] + 2
##' S2[1:g1, labels==1] <- S2[1:g1, labels==1] - 2
##' S2[g1 + 1:g2, labels==2] <- S2[g1 + 1:g2, labels==2] - 2
##' S2[g1 + 1:g2, labels==3] <- S2[g1 + 1:g2, labels==3] + 2
##' 
##' heatmap(S2, Rowv=NA, Colv=NA)
##' 
##' S = list(t(S1),t(S2))
##' groups <- Map('c',1:g1,g1 + 1:g2)
##' 
##' res <- MISKmeans(d = S, K = 3, gamma = 0.4, group = groups)
##' 
##' names(res)
##' 
##' order_S1 <- order(res$Cs[[1]])
##' order_S2 <- order(res$Cs[[2]])
##' 
##' S1_sub <- S1[res$ws!=0, order_S1]
##' S2_sub <- S2[res$ws!=0, order_S2]
##' 
##' col_S1 <- palette()[res$Cs[[1]][order_S1]]
##' col_S2 <- palette()[res$Cs[[2]][order_S2]]
##' 
##' heatmap(S1_sub, Rowv=NA, Colv=NA, ColSideColors=col_S1)
##' heatmap(S2_sub, Rowv=NA, Colv=NA, ColSideColors=col_S2)
##' 
##' 
##' 
MISKmeans <- function(d, K = NULL, gamma = NULL, lambda = 0.5, alpha = 0.5, group = NULL, nstart = 20, 
    wsPre = NULL, iniWbound = 20, penaltyInfo = NULL, silent = FALSE, maxiter = 20, sampleSizeAdjust = FALSE) {
    
    ## check input
    if (length(d) < 2) {
        stop("length of x must be greater or equal to 2.")
    }
    if (!all(sapply(d, ncol) == ncol(d[[1]]))) {
        stop("all studies must have equal number of genes (ncol)")
    }
    if (is.null(K)) {
        stop("must specify number of clusters K")
    }
    if (K < 2) {
        stop("number of clusters K must be greater than 2")
    }
    if (!is.null(penaltyInfo)) {
        if (!(length(gamma) == length(penaltyInfo))) {
            stop("gamma and penaltyInfo must have the same length.")
        }
    }
    
    
    ## obtain basic information
    numStudies <- length(d)
    J <- ncol(d[[1]])
    G0 <- length(group)
    tss.x <- list()
    
    for (i in 1:numStudies) {
        ## total sum of square for each study
        tss.x[[i]] <- apply(scale(d[[i]], center = TRUE, scale = FALSE)^2, 2, sum)
    }
    
	if(!silent){
	    cat("MISKmeans: Perform MetaSparseKmeans to initialize the result\n")	
	}
    set.seed(32608)
    mskm <- MetaSparseKmeans(d, K = K, wbounds = iniWbound, wsPre = wsPre, sampleSizeAdjust = sampleSizeAdjust, 
        silence = silent)
    # mskm1 <- MetaSpaKmeans(d, K = K, wbounds = 25, wsPre = wsPre, sampleSizeAdjust =
    # sampleSizeAdjust) mskm1 <- MetaSparseKmeans::MetaSparseKmeans(d, K = K, wbounds = 25, wsPre
    # = wsPre, sampleSizeAdjust = sampleSizeAdjust) Map(adjustedRandIndex, mskm$Cs, label)
    # Map(adjustedRandIndex, mskm1$Cs, label) Map(adjustedRandIndex, res$Cs, label) sum(mskm$ws
    # != 0)
    wsPre <- mskm$ws
    Cs <- mskm$Cs
    
    ## iteratively update CS, WS
    out <- replicate(length(gamma), list())
    for (i in 1:length(gamma)) {
        agamma <- gamma[i]
        if (is.null(penaltyInfo)) {
			if(!silent){
	            cat("Initilizaing results using alpha = 1","\n","This step will estimate the group weight penalty using unbiased feature selection principle.","\n")				
			}
            groupInfoIni <- prepareGroup(group, J, G0, agamma, 1, wsPre)
            ADMMobjectIni <- updateMISKmeans(d, K, groupInfoIni, Cs, wsPre, tss.x, lambda, sampleSizeAdjust = sampleSizeAdjust, silent=silent)
			if(!silent){
				cat("Performing MISKmeans.\n")
			}
			
			if(all(ADMMobjectIni$ws==0)){
				warning("gamma is ", agamma, ", too large such that all ws are zeor. Please consider smaller gamma.")
				out[[i]] <- ADMMobjectIni
				next
			}
			            
            groupInfo <- prepareGroup(group, J, G0, agamma, alpha, ADMMobjectIni$ws)
            ADMMobject <- updateMISKmeans(d, K, groupInfo, ADMMobjectIni$Cs, ADMMobjectIni$ws, 
                tss.x, lambda, sampleSizeAdjust = sampleSizeAdjust, silent=silent)
            # Map(adjustedRandIndex, ADMMobject$Cs, label)
        } else {
			if(!silent){
	            cat("MISKmeans: using defined groups\n")
			}
            groupInfo <- penaltyInfo[[i]]
            ADMMobject <- updateMISKmeans(d, K, groupInfo, Cs, wsPre, tss.x, lambda, sampleSizeAdjust = sampleSizeAdjust, silent=silent)
        }
        out[[i]] <- ADMMobject
    }
    
    if (length(gamma) == 1) 
        out <- out[[i]]
    class(out) <- "MISKmeans"
	if(!silent){
	    cat("MISKmeans: done\n")
	}
    return(out)
}




if(F){
	## for debug purpose
	#d=ad
	K = 3
	agamma = 1.07
	gamma=gammas
	lambda = 0.5
	alpha = 0.5
	group=groupInfo
	nstart = 20
	wsPre = NULL
	iniWbound = 20
	penaltyInfo = NULL
	silent = FALSE
	maxiter = 20
	sampleSizeAdjust = FALSE
}






