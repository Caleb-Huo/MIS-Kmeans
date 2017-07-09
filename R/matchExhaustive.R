matchExhaustive <- function(x, Cs, ws, silence = FALSE) {
    numS = length(Cs)
    corPre = prepareMCC(x, Cs)
    lenCs = vector("numeric", numS)
    uniCs = vector("list", numS)
    
    for (s in 1:numS) {
        uniCs[[s]] = unique(Cs[[s]])
        lenCs[s] = length(uniCs[[s]])
    }
    resCs = lapply(uniCs, sort)  ## result Cluster label, start form 12345, then exhausive search
    
    ## how to calculate energy, resCs specify matching rule
    highEng = eng_cor_total(corPre, reduCs = resCs, ws = ws)
    
    ### exhaustive search start here
    stdMatch = 1:lenCs[1]
    combRule = lapply(lenCs, function(x) as.matrix(combn(stdMatch, x)))
    permRule = lapply(lenCs, permn)
    combFlag = rep(1, numS)
    permFlag = rep(1, numS)
    combEndFlag = sapply(combRule, ncol)
    permEndFlag = sapply(permRule, listLength)
    tmpCombFlag = combFlag
    tmpPermFlag = permFlag
    tmpReduCs = resCs
    trace = highEng
    while (tmpPermFlag[1] == 1) {
        combFlag = tmpCombFlag
        permFlag = tmpPermFlag
        tmpEng = eng_cor_total(corPre = corPre, reduCs = tmpReduCs, ws = ws)
        #print(tmpEng)
        
        if (tmpEng > highEng) {
            highEng = tmpEng
            resCs = tmpReduCs
        }
        if (!silence) 
            # print(highEng)
        trace = c(trace, highEng)
        tmpPermFlag[[length(Cs)]] = tmpPermFlag[[length(Cs)]] + 1
        for (s in length(Cs):2) {
            if (tmpPermFlag[[s]] > permEndFlag[[s]]) {
                tmpPermFlag[[s]] = 1
                tmpCombFlag[[s]] = tmpCombFlag[[s]] + 1
                if (tmpCombFlag[[s]] > combEndFlag[[s]]) {
                  tmpCombFlag[[s]] = 1
                  tmpPermFlag[[s - 1]] = tmpPermFlag[[s - 1]] + 1
                }
            }
            tmpOrder = combRule[[s]][, tmpCombFlag[[s]]][permRule[[s]][[tmpPermFlag[s]]]]
            tmpReduCs[[s]] = reorderLabel(resCs[[s]], tmpOrder)
        }
    }
    
    #################### stop here
    perEng = eng_cor_per(corPre = corPre, reduCs = resCs)
    resumeCs = Cs
    for (s in 1:numS) {
        resumeCs[[s]] = reorderLabel(Cs[[s]], resCs[[s]])
    }
    ## return resumed matching Cs, high energy, energy per gene, trace
    return(list(matchCs = resumeCs, highEng = highEng, perEng = perEng, trace = trace))
}
