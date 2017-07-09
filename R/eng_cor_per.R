eng_cor_per <- function(corPre, reduCs) {
    p = ncol(corPre[[1]]$MeanInd)
    perEng = vector("numeric", p)
    numS = length(corPre)
    reduXComb = as.matrix(combn(numS, 2))
    
    for (i in 1:ncol(reduXComb)) {
        index1 = reduXComb[, i][1]
        index2 = reduXComb[, i][2]
        perEng = perEng + eng_MCC_pair(corPre[[index1]], corPre[[index2]], reduCs[[index1]], reduCs[[index2]])
    }
    ## scale to comparable with bcss/tss
    return(perEng/choose(length(corPre), 2))
}
