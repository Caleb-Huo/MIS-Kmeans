UpdateCs <- function(x, K, ws, Cs, nstart = nstart, tss.x) {
    newCs <- list()
    for (i in 1:length(x)) {
        x[[i]] <- x[[i]][, ws != 0]
        z <- sweep(x[[i]], 2, sqrt((ws/tss.x[[i]])[ws != 0]), "*")
        nrowz <- nrow(z)
        mus <- NULL
        if (!is.null(Cs[[i]])) {
            for (k in unique(Cs[[i]])) {
                if (sum(Cs[[i]] == k) > 1) 
                  mus <- rbind(mus, apply(z[Cs[[i]] == k, ], 2, mean))
                if (sum(Cs[[i]] == k) == 1) 
                  mus <- rbind(mus, z[Cs[[i]] == k, ])
            }
        }
        if (is.null(mus)) {
            km <- kmeans(z, centers = K, nstart = nstart)
        } else {
            distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz + 1):(nrowz + K)]
            nearest <- apply(distmat, 1, which.min)
            if (length(unique(nearest)) == K) {
                km <- kmeans(z, centers = mus)
            } else {
                km <- kmeans(z, centers = K, nstart = nstart)
            }
        }
        newCs[[i]] <- km$cluster
    }
    return(newCs)
}
