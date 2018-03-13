updateMISKmeans <- function(d, K, groupInfo, Cs, ws, tss.x, lambda, sampleSizeAdjust = FALSE, silent = FALSE, 
    maxiter = 20) {
    J <- ncol(d[[1]])
    ws.old <- rnorm(J)
    niter <- 0
    currentY <- NULL
    
    while ((sum(abs(ws - ws.old))/sum(abs(ws.old))) > 1e-04 && niter < maxiter) {
        if (!silent) 
            cat("Iteration", niter, ":\n", fill = FALSE)
        niter <- niter + 1
        ws.old <- ws
        if (sum(ws != 0) < 1) {
            wsPre <- ws
            objective <- 0
            obj0 <- 0
            break
        }
        if (!silent) 
            cat("Updating CS...\n", fill = FALSE)
        if (niter > 1) 
            Cs <- UpdateCs(d, K, ws, Cs, tss.x)  # if niter=1, no need to update!!
        ## UpdateCs(d, K, ws, Cs, tss.x)
        if (!silent) 
            cat("Updating WS...\n", fill = FALSE)
        if (is.null(groupInfo)) {
            fmatch = patternMatch(d, Cs, ws, silence = silence)
            ratio = GetRatio(d, Cs, tss.x, sampleSizeAdjust = sampleSizeAdjust)
            aa <- ratio + lambda * (fmatch$perEng + 1)/2
            
            ws <- aa/sqrt(sum(aa^2))
            objective <- -sum(ws * aa)
            obj0 <- -sum(ws * aa)
            # print(objective)
            
            Cs_match <- fmatch$Cs
            
        } else {
            ADMMobject <- UpdateWsADMM_m(d, Cs, ws, currentY = currentY, groupInfo, tss.x, lambda, sampleSizeAdjust = sampleSizeAdjust)
            ws <- ADMMobject$z
            # print(sum(ws != 0))
            currentY <- ADMMobject$currentY
            # print(ADMMobject$objective)
            
            objective = ADMMobject$objective
            obj0 <- ADMMobject$obj0
            
            Cs_match <- ADMMobject$Cs
        }
    }
    
    res <- list(ws = ws, Cs = Cs_match, obj0 = obj0, objective = objective, groupInfo = groupInfo)
    return(res)
}


UpdateWsADMM_m <- function(d, Cs, ws, currentY = NULL, groupInfo, tss.x, lambda, sampleSizeAdjust = FALSE) {
    
    fmatch = patternMatch(d, Cs, ws, silence = silence)
    ratio = GetRatio(d, Cs, tss.x, sampleSizeAdjust = sampleSizeAdjust)
    aa <- ratio + lambda * (fmatch$perEng + 1)/2
    
    J <- groupInfo$J
    L <- groupInfo$L
    G <- groupInfo$G
    groupLevel <- groupInfo$groupLevel
    genePos <- groupInfo$genePos - 1
    coef <- groupInfo$coef
    
    
    stopifnot(L == length(groupLevel))
    stopifnot(L == length(genePos))
    stopifnot(L == length(coef))
    stopifnot(1:L == order(groupLevel))
    stopifnot(max(groupLevel) == G)
    stopifnot(max(genePos) == J - 1)
    
    if (is.null(currentY)) 
        currentY <- numeric(L)
    
    x <- numeric(L)
    z <- ws
    
    ADMMobj <- .C("ADMM_updatew_R", x = as.double(x), currentY = as.double(currentY), z = as.double(z), r = as.double(aa), 
        objective = as.double(0), groupLevel = as.integer(groupLevel), genePos = as.integer(genePos), coef = as.double(coef), 
        J = as.integer(J), G = as.integer(G), L = as.integer(L))
    
    ADMMobj$x <- NULL
    ADMMobj$r <- NULL
    
    ADMMobj$groupLevel <- NULL
    ADMMobj$genePos <- NULL
    ADMMobj$coef <- NULL
    ADMMobj$J <- NULL
    ADMMobj$G <- NULL
    ADMMobj$L <- NULL
    
    ADMMobj$obj0 <- -sum(ws * aa)
    ADMMobj$Cs = fmatch$matchCs
    
    return(ADMMobj)
}


## prepare group information here module should be G-lists. Each element g contain feature indexes whose
## domain is P.  L is the expanded length of non-zero element in J by G0 matrix.  groupLevel (L):
## 1,1,1,1,2,2,2,3,3,3,... increasing, same number indicate same group.  genePos (L): position.  coef (L):
## coef for the expanded features.  z (J): is the feature weight.  x (J): primal variable.  y (J): dual
## variable.  ws: feature weight of previous iteration

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
            # cat(agroupPenalty) cat(' ')
            coef[curPos:endPos] <- preCoef * sqrt(a_inv_groupFeatureCounts) * sqrt(agroupPenalty)
            curPos <- curPos + alen
        }
        # cat('\n')
        
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
            # cat(agroupPenalty) cat(' ')
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
    
    groupInfo <- list(groupLevel = groupLevel, genePos = genePos, coef = coef, L = L, G = G0 + J0, J = J, alpha = alpha, 
        gamma = gamma)
    return(groupInfo)
}



