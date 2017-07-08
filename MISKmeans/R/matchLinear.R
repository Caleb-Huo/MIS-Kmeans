matchLinear <- function(x, Cs, ws, silence = FALSE) {
    numS = length(Cs)
    corPre = prepareMCC(x, Cs)
    lenCs = vector("numeric", numS)
    uniCs = vector("list", numS)
    
    studySize <- sapply(Cs, length)
    studySizeRank <- rank(studySize, ties.method = "random")
    
    for (s in 1:numS) {
        uniCs[[s]] = unique(Cs[[s]])
        lenCs[s] = length(uniCs[[s]])
    }
    
    ## initialize
    resCs <- uniCs
    
    
    for (i in length(studySizeRank):1) {
        if (i == length(studySizeRank)) {
            thisIndex <- which(studySizeRank == i)
            currentIndexes <- thisIndex
            resCs[[thisIndex]] = sort(resCs[[thisIndex]])
        } else {
            thisIndex <- which(studySizeRank == i)
            currentIndexes <- c(currentIndexes, thisIndex)
            thiscorPre <- corPre[currentIndexes]
            thisresCs <- resCs[currentIndexes]
            tmpThisResCs <- thisresCs
            highEng <- eng_cor_total(thiscorPre, reduCs = thisresCs, ws = ws)
            
            permRule = permn(lenCs[thisIndex])
            innerThisIndex <- length(currentIndexes)
            
            for (j in 1:length(permRule)) {
                tmpThisResCs[[innerThisIndex]] = permRule[[j]]
                thisEnergy <- eng_cor_total(thiscorPre, reduCs = tmpThisResCs, ws = ws)
                if (thisEnergy > highEng) {
                  highEng <- thisEnergy
                  thisresCs <- tmpThisResCs
                  resCs[[thisIndex]] <- permRule[[j]]
                }
            }
        }
    }
    
    #################### stop here
    perEng = eng_cor_per(corPre = corPre, reduCs = resCs, ws = ws)
    resumeCs = Cs
    for (s in 1:numS) {
        resumeCs[[s]] = reorderLabel(Cs[[s]], resCs[[s]])
    }
    ## return resumed matching Cs, high energy, energy per gene, trace
    return(list(matchCs = resumeCs, highEng = highEng, perEng = perEng, trace = NULL))
    
}
