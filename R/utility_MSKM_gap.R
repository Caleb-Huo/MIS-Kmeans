permuteX <- function(x) {
    permx <- list()
    for (i in 1:length(x)) {
        permx[[i]] <- matrix(NA, nrow = nrow(x[[i]]), ncol = ncol(x[[i]]))
        for (j in 1:ncol(x[[i]])) permx[[i]][, j] <- sample(x[[i]][, j])
    }
    return(permx)
}
