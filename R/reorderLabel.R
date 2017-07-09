reorderLabel <- function(alabel, aorder) {
    uniLab = unique(alabel)
    resLabel = NULL
    for (i in 1:length(uniLab)) {
        resLabel[uniLab[i] == alabel] = aorder[i]
    }
    return(resLabel)
}
