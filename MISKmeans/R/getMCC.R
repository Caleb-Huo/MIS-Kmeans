getMCC <- function(S1, S2, label1, label2) {
    S1 = as.matrix(S1)
    S2 = as.matrix(S2)
    muTot1 = colMeans(S1)
    muTot2 = colMeans(S2)
    p = ncol(S1)
    
    uniLabel = unique(label1)
    K = length(uniLabel)
    
    muInd1 = matrix(NA, nrow = K, ncol = p)
    muInd2 = matrix(NA, nrow = K, ncol = p)
    varInd1 = matrix(NA, nrow = K, ncol = p)
    varInd2 = matrix(NA, nrow = K, ncol = p)
    
    for (i in 1:K) {
        muInd1[i, ] = colMeans(S1[i == label1, ])
        muInd2[i, ] = colMeans(S2[i == label2, ])
        varInd1[i, ] = apply(S1[i == label1, ], 2, var)
        varInd2[i, ] = apply(S2[i == label2, ], 2, var)
    }
    
    MCC_num = colMeans(muInd1 * muInd2) - colMeans(muInd1) * colMeans(muInd2)
    den1 = colMeans(varInd1) + colMeans((muInd1 - matrix(rep(muTot1, each = K), nrow = K))^2)
    den2 = colMeans(varInd2) + colMeans((muInd2 - matrix(rep(muTot2, each = K), nrow = K))^2)
    MCC_den = sqrt(den1 * den2)
    MCC = MCC_num/MCC_den
    return(MCC)
}
