UpdateCs1 <- function(x, K, ws, Cs, tss.x, nstart = nstart) {
    x <- x[, ws != 0]
    z <- sweep(x, 2, sqrt((ws/tss.x)[ws != 0]), "*")
    nrowz <- nrow(z)
    mus <- NULL
    if (!is.null(Cs)) {
        for (k in unique(Cs)) {
            if (sum(Cs == k) > 1) 
                mus <- rbind(mus, apply(z[Cs == k, ], 2, mean))
            if (sum(Cs == k) == 1) 
                mus <- rbind(mus, z[Cs == k, ])
        }
    }
    if (is.null(mus)) {
        km <- kmeans(z, centers = K, nstart = nstart)
    } else {
        distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz + 1):(nrowz + K)]
        nearest <- apply(distmat, 1, which.min)
        if (length(unique(nearest)) == K) {
            km <- kmeans(z, centers = mus)
        } else {
            km <- kmeans(z, centers = K, nstart = nstart)
        }
    }
    newCs <- km$cluster
    return(newCs)
}


GetWCSS <- function(x, Cs, ws = NULL) {
    wcss.perfeature <- numeric(ncol(x))
    for (k in unique(Cs)) {
        whichers <- (Cs == k)
        if (sum(whichers) > 1) 
            wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers, ], center = TRUE, scale = FALSE)^2, 
                2, sum)
    }
    tss.perfeature <- apply(scale(x, center = TRUE, scale = FALSE)^2, 2, sum)
    bcss.perfeature <- tss.perfeature - wcss.perfeature
    r <- bcss.perfeature/tss.perfeature
    
    if (!is.null(ws)) 
        return(list(wcss.perfeature = wcss.perfeature, wcss = sum(wcss.perfeature), wcss.ws = sum(wcss.perfeature * 
            ws), bcss.perfeature = bcss.perfeature, r = r))
    if (is.null(ws)) 
        return(list(wcss.perfeature = wcss.perfeature, wcss = sum(wcss.perfeature), bcss.perfeature = bcss.perfeature, 
            r = r))
}



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
            print(objective)
        } else {
            ADMMobject <- UpdateWsADMM_m(d, Cs, ws, currentY = currentY, groupInfo, tss.x, lambda, sampleSizeAdjust = sampleSizeAdjust)
            ws <- ADMMobject$z
            # print(sum(ws != 0))
            currentY <- ADMMobject$currentY
            print(ADMMobject$objective)
            
            objective = ADMMobject$objective
            obj0 <- ADMMobject$obj0
        }
    }
    
    res <- list(ws = ws, Cs = ADMMobject$Cs, obj0 = obj0, objective = objective, groupInfo = groupInfo)
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
    
    ADMMobj <- .C("ADMM_updatew_R", x = as.double(x), currentY = as.double(currentY), z = as.double(z), 
        r = as.double(aa), objective = as.double(0), groupLevel = as.integer(groupLevel), genePos = as.integer(genePos), 
        coef = as.double(coef), J = as.integer(J), G = as.integer(G), L = as.integer(L))
    
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


