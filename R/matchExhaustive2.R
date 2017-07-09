matchExhaustive2 <- function(x, Cs, ws, silence = FALSE) {
    numS = length(Cs)
    
    ws2 <- ws[ws != 0]
    
    corPre = prepareMCC(x, Cs)
    corPre2 <- corPre
	for(s in 1:length(corPre2)){
		corPre2[[s]]$MeanInd <- corPre2[[s]]$MeanInd[,ws2]
		corPre2[[s]]$denInd <- corPre2[[s]]$denInd[ws2]
	}
		
    S <- length(corPre)
    G <- length(ws)
    
    CombS <- as.data.frame(combn(S, 2))
    encodeS <- sapply(CombS, encode, S)
    
    K <- length(unique(Cs[[1]]))
    
    permK <- permn(K)
    encodeK <- sapply(permK, encode, K)
    
    
    energyS <- replicate(length(CombS), list())
    for (i in seq_along(CombS)) {
        index_a <- CombS[[i]][1]
        index_b <- CombS[[i]][2]
        
        corPre_a <- corPre[[index_a]]
        corPre_b <- corPre[[index_b]]
        
        energyK <- replicate(length(permK), list())
        for (j in seq_along(permK)) {
            bRank <- permK[[j]]
            energyK[[j]] <- sum(eng_MCC_pair(corPre_a, corPre_b, 1:K, bRank) * ws2)
        }
        ahashK <- hash(encodeK, energyK)
        energyS[[i]] <- ahashK
    }
    if (length(energyS) == 1) {
        energyS <- energyS[[1]]
    }
    
    hashS <- hash(encodeS, energyS)
    
    resCs <- as.data.frame(replicate(S, 1:K))
    
    ### exhaustive search start here
    lenCs <- rep(K, S)
    permRule = lapply(lenCs, permn)
    permFlag = rep(1, S)
    permEndFlag = sapply(permRule, listLength)
    
    tmpCs = resCs
    
    iniEnergy <- numeric(S)
    for (s in 1:S) {
        if (s > 1) {
            interEnergy <- 0
            for (ps in 1:(s - 1)) {
                ps_order <- order(tmpCs[[ps]])
                Cs_ps <- tmpCs[[ps]][ps_order]
                aSencode <- as.character(encode(c(ps, s), S))
                
                Cs_s <- tmpCs[[s]][ps_order]
                aKencode <- as.character(encode(Cs_s, K))
                
                interEnergy <- interEnergy + hashS[[aSencode]][[aKencode]]
                
            }
            iniEnergy[s] <- iniEnergy[s - 1] + interEnergy
        }
    }
    
    highEng <- iniEnergy[S]
    
    while (permFlag[1] == 1) {
        
        tmpCs[[S]] <- permRule[[S]][[permFlag[[S]]]]
        
        interEnergy <- 0
        for (ps in 1:(S - 1)) {
            ps_order <- order(tmpCs[[ps]])
            Cs_ps <- tmpCs[[ps]][ps_order]
            aSencode <- as.character(encode(c(ps, S), S))
            
            Cs_s <- tmpCs[[S]][ps_order]
            aKencode <- as.character(encode(Cs_s, K))
            interEnergy <- interEnergy + hashS[[aSencode]][[aKencode]]
        }
        
        tmpEng <- iniEnergy[S - 1] + interEnergy
        #print(tmpEng/choose(S, 2))
        if (tmpEng > highEng) {
            highEng = tmpEng
            resCs = tmpCs
        }
        permFlag[[S]] = permFlag[[S]] + 1
        
        updateFlag <- rep(0, S)
        for (s in S:2) {
            if (permFlag[[s]] > permEndFlag[[s]]) {
                permFlag[[s]] = 1
                permFlag[[s - 1]] = permFlag[[s - 1]] + 1
                updateFlag[s] <- 1
                updateFlag[s - 1] <- 1
                
            }
        }
        
        for (s in S:2) {
            if (updateFlag[s] == 1) {
                tmpCs[[s]] <- permRule[[s]][[permFlag[[s]]]]
                
                interEnergy <- 0
                for (ps in 1:(s - 1)) {
                  ps_order <- order(tmpCs[[ps]])
                  Cs_ps <- tmpCs[[ps]][ps_order]
                  aSencode <- as.character(encode(c(ps, s), S))
                  
                  Cs_s <- tmpCs[[s]][ps_order]
                  aKencode <- as.character(encode(Cs_s, K))
                  interEnergy <- interEnergy + hashS[[aSencode]][[aKencode]]
                }
                
                iniEnergy[s:S] <- iniEnergy[s:S] + interEnergy - iniEnergy[s] + iniEnergy[s - 1]
            }
        }
    }
    
    
    #################### stop here
    perEng = eng_cor_per(corPre = corPre, reduCs = resCs)
    resumeCs = Cs
    for (s in 1:numS) {
        resumeCs[[s]] = reorderLabel(Cs[[s]], resCs[[s]])
    }
    ## return resumed matching Cs, high energy, energy per gene, trace
    return(list(matchCs = resumeCs, highEng = highEng, perEng = perEng))
}
