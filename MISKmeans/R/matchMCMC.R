matchMCMC <- function(x, Cs, ws, silence = FALSE) {
    
    matchLinearRes <- matchLinear(x, Cs, ws, silence = silence)
    thisCs = matchLinearRes$matchCs  ## result Cluster label, start form 12345, then exhausive search
    highEng = matchLinearRes$highEng
    
    numS = length(thisCs)
    corPre = prepareMCC(x, thisCs)
    lenCs = vector("numeric", numS)
    uniCs = vector("list", numS)
    
    for (s in 1:numS) {
        uniCs[[s]] = unique(thisCs[[s]])
        lenCs[s] = length(uniCs[[s]])
    }
    resCs <- uniCs
    ############# start here
    
    sim_num = 300
    step_max = 1e+05
    count = 0
    highestCs = tmpReduCs = resCs
    highestEng = T_schedule = highEng
    trace = highestEng
    
    stdMatch = 1:lenCs[1]
    
    while (count < step_max) {
        accept = 0
        # print(paste('temperature:',T_schedule))
        for (t in 1:sim_num) {
            count = count + 1
            study_index = sample(1:length(x), 1)
            element_index = sample(stdMatch, 2)
            tmpReduCs[[study_index]] = rearrange(resCs[[study_index]], element_index)
            tmpEng = eng_cor_total(corPre, reduCs = tmpReduCs, ws = ws)
            if (tmpEng > highEng) {
                accept = accept + 1
                highEng = tmpEng
                resCs = tmpReduCs
            } else {
                prob = exp(-(highEng - tmpEng)/T_schedule)
                if (runif(1) < prob) {
                  accept = accept + 1
                  highEng = tmpEng
                  resCs = tmpReduCs
                }
            }
            trace = c(trace, highEng)
            if (highestEng < highEng) {
                highestEng = highEng
                highestCs = resCs
            }
        }
        
        
        acc_ratio = accept/sim_num
        # print(acc_ratio)
        if (acc_ratio < 0.1) {
            #### stopping rule
            break
        } else {
            if (acc_ratio > 0.7) {
                T_schedule = T_schedule * 0.7
            } else {
                T_schedule = T_schedule * 0.9
            }
        }
    }
    
    highEng = highestEng
    resCs = highestCs
    #################### stop here
    
    perEng = eng_cor_per(corPre = corPre, reduCs = resCs, ws = ws)
    resumeCs = Cs
    for (s in 1:numS) {
        resumeCs[[s]] = reorderLabel(Cs[[s]], resCs[[s]])
    }
    ## return resumed matching Cs, high energy, energy per gene, trace
    return(list(matchCs = resumeCs, highEng = highEng, perEng = perEng, trace = trace))
}
