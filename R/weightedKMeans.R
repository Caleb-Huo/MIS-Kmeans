weightedKMeans <- function(x, K, ws, tss.x = NULL) {
    if (is.null(tss.x)) {
        tss.x <- apply(scale(t(x), center = TRUE, scale = FALSE)^2, 2, sum)
    }
    
    commonNonZeroNames = intersect(rownames(x), names(ws)[ws != 0])
    x <- x[commonNonZeroNames, ]
    tss.x = tss.x[commonNonZeroNames]
    z <- sweep(x, 1, sqrt(ws[commonNonZeroNames]/tss.x), "*")
    km <- kmeans(t(z), centers = K, nstart = 50)
    newCs <- km$cluster
    return(newCs)
}
