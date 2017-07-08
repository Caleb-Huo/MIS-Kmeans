eng_future_total <- function(corPreStar, reduCsStar, wsStar, corPreNew, reduCsNew) {
    
    perEng = vector("numeric", length(wsStar))
    
    numS = length(corPreStar)
    
    for (i in 1:numS) {
        perEng = perEng + eng_MCC_pair(corPreStar[[i]], corPreNew[[1]], reduCsStar[[i]], reduCsNew)
    }
    ## scale to comparable with bcss/tss
    perEng = perEng/numS
    return(sum(perEng * wsStar))
}
