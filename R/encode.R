encode <- function(avec, n) {
    res <- 0
    for (i in seq_along(avec)) {
        res <- res + avec[i] * n^(i - 1)
    }
    res
}

