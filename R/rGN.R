rGN <- function(ax, alabel) {
    level = unique(alabel)
    len = length(level)
    M = numeric(len)
    SD = numeric(len)
    for (i in 1:length(level)) {
        M[i] = mean(ax[level[i] == alabel])
        SD[i] = sd(ax[level[i] == alabel])
    }
    a_rGN = mean(M)
    b_rGN = mean(SD^2) + mean((M - a_rGN)^2)
    res = ((ax - a_rGN)/sqrt(b_rGN))
    return(res)
}
