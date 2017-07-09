GetRatio <- function(x, Cs, tss.x, sampleSizeAdjust = FALSE) {
    ratio.bcss_tss <- numeric(ncol(x[[1]]))
    for (i in 1:length(x)) {
        wcss.perfeature <- numeric(ncol(x[[i]]))
        tss.perfeature <- tss.x[[i]]
        for (k in unique(Cs[[i]])) {
            whichers <- (Cs[[i]] == k)
            if (sum(whichers) > 1) 
                wcss.perfeature <- wcss.perfeature + apply(scale(x[[i]][whichers, ], center = TRUE, scale = FALSE)^2, 
                  2, sum)
        }
        aratio <- numeric(ncol(x[[1]]))
        bcss.perfeature = tss.perfeature - wcss.perfeature
        tssNonZeroIndex <- tss.perfeature != 0
        aratio[tssNonZeroIndex] <- bcss.perfeature[tssNonZeroIndex]/tss.perfeature[tssNonZeroIndex]
        
        if (sampleSizeAdjust) {
            thisNsamples <- sapply(x, function(x) nrow(x))
            ratio.bcss_tss = ratio.bcss_tss + aratio * thisNsamples[i]/sum(thisNsamples)
        } else {
            ratio.bcss_tss = ratio.bcss_tss + aratio/length(x)
        }
    }
    ratio.bcss_tss
}
