prepareMCC <- function(x, Cs) {
    numS = length(Cs)
    p = ncol(x[[1]])
    res = list()
    for (i in 1:numS) {
        reduCs = unique(Cs[[i]])
        m = length(reduCs)
        reduMeanTot = colMeans(x[[i]])
        MeanInd = matrix(NA, nrow = m, ncol = p)
        VarInd = matrix(NA, nrow = m, ncol = p)
        MeanTot = colMeans(x[[i]])
        for (j in 1:m) {
            tmpVar = x[[i]][reduCs[j] == Cs[[i]], ]
            if (is.null(dim(tmpVar))) {
                MeanInd[j, ] = tmpVar
                VarInd[j, ] = 0
            } else {
                MeanInd[j, ] = colMeans(tmpVar)
                VarInd[j, ] = apply(tmpVar, 2, var)
            }
            zeroIndex = VarInd[j, ] == 0
            # VarInd[j,zeroIndex] = rep(0.001,sum(zeroIndex))
            VarInd[j, ][zeroIndex] = 1e-04
        }
        denInd = colMeans(VarInd) + colMeans((MeanInd - matrix(rep(MeanTot, each = m), nrow = m))^2)
        res[[i]] = list(MeanInd = MeanInd, denInd = denInd)
    }
    return(res)
}
