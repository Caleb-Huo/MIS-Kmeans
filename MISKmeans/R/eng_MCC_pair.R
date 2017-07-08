eng_MCC_pair <- function(corS1, corS2, reduCs1, reduCs2) {
    K1 = length(reduCs1)
    K2 = length(reduCs2)
    p = ncol(corS1$MeanInd)
    m = K1
    EXY = numeric(p)
    EX = numeric(p)
    EY = numeric(p)
    
    for (i in 1:m) {
        EXY = EXY + corS1$MeanInd[reduCs1 == i, ] * corS2$MeanInd[reduCs2 == i, ]
        EX = EX + corS1$MeanInd[reduCs1 == i, ]
        EY = EY + corS2$MeanInd[reduCs2 == i, ]
    }
    MCC_num = EXY/m - EX * EY/m^2
    MCC_den = sqrt(corS1$denInd * corS2$denInd)
    MCC_pair_per = MCC_num/MCC_den
    return(MCC_pair_per)
}
