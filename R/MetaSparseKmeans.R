MetaSparseKmeans <- function(x, K = NULL, wbounds = NULL, nstart = 20, ntrial = 1, maxiter = 20, lambda = 1/2, 
    sampleSizeAdjust = FALSE, wsPre = NULL, silence = FALSE) {
    
    ## check input
    if (length(x) < 2) {
        stop("length of x must be greater or equal to 2.")
    }
    if (!all(sapply(x, ncol) == ncol(x[[1]]))) {
        stop("all studies must have equal number of genes (ncol)")
    }
    if (is.null(K)) {
        stop("must specify number of clusters K")
    }
    if (K < 2) {
        stop("number of clusters K must be greater than 2")
    }
    
    
    ## get some basic information
    numStudies = length(x)
    
    ## create result
    bestOut <- vector("list", length(wbounds))
    out <- vector("list", length(wbounds))
    Cs0 <- list()
    tss.x <- list()
    
    for (i in 1:numStudies) {
        ## total sum of square for each study
        tss.x[[i]] <- apply(scale(x[[i]], center = TRUE, scale = FALSE)^2, 2, sum)
    }
    
    for (atrail in 1:ntrial) {
        # initialize initialize cluster by KMeans initialize w
        if (is.null(wsPre)) {
            wsPre <- numeric(ncol(x[[1]]))
            for (i in 1:numStudies) {
                asparcl <- KMeansSparseCluster(x[[i]], K = K, wbounds = wbounds[1], silent = silence)[[1]]
                Cs0[[i]] <- asparcl$Cs
                wsPre <- wsPre + asparcl$ws/numStudies
            }
        } else {
            if (length(wsPre) != ncol(x[[1]])) 
                stop("length of wsPre differs from number of genes")
            if (is.null(names(wsPre))) 
                stop("there is no name for wsPre")
            if (any(names(wsPre) != colnames(x[[1]]))) 
                stop("name of wsPre differs from gene name")
            for (i in 1:numStudies) Cs0[[i]] <- weightedKMeans(x = t(x[[i]]), K = K, ws = wsPre, tss.x = tss.x[[i]])
        }
        # 
        for (w in 1:length(wbounds)) {
            awbound = wbounds[w]
            
            ws <- wsPre
            
            ws.old <- rnorm(ncol(x[[1]]))
            store.ratio <- NULL
            niter <- 0
            trace <- list()
            Cs = Cs0
            
            ## iterate until converge or maxiteration
            while ((sum(abs(ws - ws.old))/sum(abs(ws.old))) > 1e-04 && niter < maxiter) {
                niter <- niter + 1
                ws.old <- ws
                if (niter > 1) 
                  Cs <- UpdateCs(x, K, ws, Cs, tss.x, nstart = nstart)  # if niter=1, no need to update!!
                
                # fmatch = patternMatch(x, Cs, ws, silence = silence)
                fmatch = patternMatch(x, Cs, ws, silence = silence)
                
                
                ratio = GetRatio(x, Cs, tss.x, sampleSizeAdjust = sampleSizeAdjust)
                ws <- UpdateWs(x, Cs, awbound, ratio, lambda * (fmatch$perEng + 1)/2)
                store.ratio <- c(store.ratio, sum(ratio * ws))
            }
            score = sum((ratio + lambda * (fmatch$perEng + 1)/2) * ws)
            names(ws) <- colnames(x[[1]])
            out[[w]] <- list(ws = ws, Cs = fmatch$matchCs, wbound = awbound, score = score)
            if (is.null(bestOut[[w]])) {
                bestOut[[w]] = out[[w]]
            } else {
                if (bestOut[[w]]$score < out[[w]]$score) 
                  bestOut[[w]] = out[[w]]
            }
        }
        
    }
    if (length(bestOut) == 1) {
        bestOut = bestOut[[1]]
    }
    class(bestOut) <- "MetaSparseKmeans"
    return(bestOut)
}
