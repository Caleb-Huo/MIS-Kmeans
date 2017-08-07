##' @useDynLib MISKmeans
##' @export
MISKmeans <- function(d, K = NULL, gamma = NULL, lambda = 0.5, alpha = 0.5, group = NULL, nstart = 20, 
    wsPre = NULL, penaltyInfo = NULL, silent = FALSE, maxiter = 20, sampleSizeAdjust = FALSE) {
    
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
    if (K < 3) {
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
    
    mskm <- MetaSparseKmeans(d, K = K, wbounds = 12, wsPre = wsPre, sampleSizeAdjust = sampleSizeAdjust)
    # Map(adjustedRandIndex, mskm$Cs, label)
    
    wsPre <- mskm$ws
    Cs <- mskm$Cs
    
    ## iteratively update CS, WS
    out <- replicate(length(gamma), list())
    for (i in 1:length(gamma)) {
        agamma <- gamma[i]
        if (is.null(penaltyInfo)) {
            cat("initilizaing results using alpha = 1\n")
            groupInfoIni <- prepareGroup(group, J, G0, agamma, 1, wsPre)
            ADMMobjectIni <- updateMISKmeans(d, K, groupInfoIni, Cs, wsPre, tss.x, lambda, sampleSizeAdjust = sampleSizeAdjust)
            cat("initilizaing groups\n")
            groupInfo <- prepareGroup(group, J, G0, agamma, alpha, ADMMobjectIni$ws)
            ADMMobject <- updateMISKmeans(d, K, groupInfo, ADMMobjectIni$Cs, ADMMobjectIni$ws, tss.x, lambda, 
                sampleSizeAdjust = sampleSizeAdjust)
            # Map(adjustedRandIndex, ADMMobject$Cs, label)
        } else {
            cat("using defined groups\n")
            groupInfo <- penaltyInfo[[i]]
            ADMMobject <- updateMISKmeans(d, K, groupInfo, Cs, wsPre, tss.x, sampleSizeAdjust = sampleSizeAdjust)
        }
        out[[i]] <- ADMMobject
    }
    
    if (length(gamma) == 1) 
        out <- out[[i]]
    class(out) <- "MISKmeans"
    return(out)
}


## prepare group information here module should be G-lists. Each element g contain feature indexes whose
## domain is P.  L is the expanded length of non-zero element in J by G0 matrix.  groupLevel (L):
## 1,1,1,1,2,2,2,3,3,3,... increasing, same number indicate same group.  genePos (L): position.  coef
## (L): coef for the expanded features.  z (J): is the feature weight.  x (J): primal variable.  y (J):
## dual variable.  ws: feature weight of previous iteration

prepareGroup <- function(group, J, G0, gamma, alpha, ws) {
    ## take care of trivial class
    if (gamma == 0) {
        return(NULL)
    }
    
    groupFeatureCounts <- numeric(J)
    for (g in 1:G0) {
        groupFeatureCounts[group[[g]]] <- groupFeatureCounts[group[[g]]] + 1
    }
    
    curPos <- 1
    preCoef <- gamma * (1 - alpha)
    
    if (alpha == 0) {
        J0logic <- groupFeatureCounts == 0
        J0 = sum(J0logic)
        
        L <- sum(groupFeatureCounts) + J0
        groupLevel <- numeric(L)
        genePos <- numeric(L)
        coef <- numeric(L)
        
        for (g in 1:G0) {
            agroup <- group[[g]]
            aws <- ws[agroup]
            alen <- length(agroup)
            endPos <- curPos + alen - 1
            groupLevel[curPos:endPos] <- g
            genePos[curPos:endPos] <- agroup
            a_inv_groupFeatureCounts <- 1/groupFeatureCounts[agroup]
            agroupPenalty <- max(sum(a_inv_groupFeatureCounts[aws != 0]), 1)
            cat(agroupPenalty)
            cat(" ")
            coef[curPos:endPos] <- preCoef * sqrt(a_inv_groupFeatureCounts) * sqrt(agroupPenalty)
            curPos <- curPos + alen
        }
        cat("\n")
        
        endPos <- curPos + J0 - 1
        groupLevel[curPos:endPos] <- (G0 + 1):(G0 + J0)
        genePos[curPos:endPos] <- (1:J)[J0logic]
        coef[curPos:endPos] <- gamma
        
    } else if (alpha == 1) {
        J0 = J
        
        L <- J
        groupLevel <- numeric(J)
        genePos <- numeric(J)
        coef <- numeric(J)
        G0 = 0
        
        endPos <- curPos + J - 1
        groupLevel[curPos:endPos] <- (G0 + 1):(G0 + J)
        genePos[curPos:endPos] <- 1:J
        coef[curPos:endPos] <- alpha * gamma
    } else {
        J0 = J
        L <- sum(groupFeatureCounts) + J
        groupLevel <- numeric(L)
        genePos <- numeric(L)
        coef <- numeric(L)
        
        for (g in 1:G0) {
            agroup <- group[[g]]
            aws <- ws[agroup]
            alen <- length(agroup)
            endPos <- curPos + alen - 1
            groupLevel[curPos:endPos] <- g
            genePos[curPos:endPos] <- agroup
            a_inv_groupFeatureCounts <- 1/groupFeatureCounts[agroup]
            agroupPenalty <- max(sum(a_inv_groupFeatureCounts[aws != 0]), 1)
            cat(agroupPenalty)
            cat(" ")
            coef[curPos:endPos] <- preCoef * sqrt(a_inv_groupFeatureCounts) * sqrt(agroupPenalty)
            curPos <- curPos + alen
        }
        cat("\n")
        
        endPos <- curPos + J - 1
        groupLevel[curPos:endPos] <- (G0 + 1):(G0 + J)
        genePos[curPos:endPos] <- 1:J
        coef[curPos:endPos] <- alpha * gamma
        coef[curPos:endPos][groupFeatureCounts == 0] <- gamma
    }
    
    groupInfo <- list(groupLevel = groupLevel, genePos = genePos, coef = coef, L = L, G = G0 + J0, J = J, 
        alpha = alpha, gamma = gamma)
    return(groupInfo)
}



