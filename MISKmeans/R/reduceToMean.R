reduceToMean <- function(S, label) {
    len = length(S)
    uniLabel = vector(mode = "list", length = len)
    resMat = vector(mode = "list", length = len)
    for (k in 1:len) {
        aS = S[[k]]
        alabel = label[[k]]
        auniLabel = unique(alabel)
        aresMat = matrix(NA, ncol = ncol(aS), nrow = length(auniLabel))
        for (i in 1:length(auniLabel)) {
            aresMat[i, ] = colMeans(aS[alabel == auniLabel[i], ])
        }
        uniLabel[[k]] = auniLabel
        resMat[[k]] = aresMat
    }
    return(list(reduMat = resMat, uniLabel = uniLabel))
}
