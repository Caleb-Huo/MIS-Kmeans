eng_cor_total <- function(corPre, reduCs, ws) {
    non0ws = ws != 0
    perEng = vector("numeric", sum(non0ws))
    
    numS = length(corPre)
    for (i in 1:numS) {
        corPre[[i]]$MeanInd = corPre[[i]]$MeanInd[, non0ws]
        corPre[[i]]$denInd = corPre[[i]]$denInd[non0ws]
    }
    reduXComb = as.matrix(combn(numS, 2))
    
    for (i in 1:ncol(reduXComb)) {
        index1 = reduXComb[, i][1]
        index2 = reduXComb[, i][2]
        perEng = perEng + eng_MCC_pair(corPre[[index1]], corPre[[index2]], reduCs[[index1]], reduCs[[index2]])
    }
    ## scale to comparable with bcss/tss
    perEng = perEng/choose(length(corPre), 2)
    return(sum(perEng * ws[non0ws]))
}
