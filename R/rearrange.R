rearrange <- function(alabel, aindex) {
    ## aindex contains two indexes to be exchanged
    alabel[which(aindex[1] == alabel)] = 0
    alabel[which(aindex[2] == alabel)] = aindex[1]
    alabel[which(0 == alabel)] = aindex[2]
    return(alabel)
}
