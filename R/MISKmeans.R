##' @export
MISKmeans <- function(d, K = NULL, gamma = NULL, lambda = 0.5, alpha = 0.5, group = NULL, 
    nstart = 20, wsPre = NULL, penaltyInfo = NULL, silent = FALSE, maxiter = 20, sampleSizeAdjust = FALSE) {
    
    if (!is.null(penaltyInfo)) {
        if (!(length(gamma) == length(penaltyInfo))) {
            stop("gamma and penaltyInfo must have the same length.")
        }
    }
    if (is.null(K)) 
        stop("Must provide either K or centers.")
    
    ## obtain basic information
    numStudies <- length(d)
    J <- ncol(d[[1]])
    G0 <- length(group)
    tss.x <- list()
    
    for (i in 1:numStudies) {
        ## total sum of square for each study
        tss.x[[i]] <- apply(scale(d[[i]], center = TRUE, scale = FALSE)^2, 2, sum)
    }
    
    mskm <- MetaSparseKmeans(d, K = K, wbounds = 12, wsPre = wsPre)
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
            ADMMobjectIni <- updateMISKmeans(d, K, groupInfoIni, Cs, wsPre, tss.x, lambda, 
                sampleSizeAdjust = sampleSizeAdjust)
            cat("initilizaing groups\n")
            groupInfo <- prepareGroup(group, J, G0, agamma, alpha, ADMMobjectIni$ws)
            ADMMobject <- updateMISKmeans(d, K, groupInfo, ADMMobjectIni$Cs, ADMMobjectIni$ws, 
                tss.x, lambda, sampleSizeAdjust = sampleSizeAdjust)
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

