BinarySearch <- function(argu, sumabs) {
    if (l2n(argu) == 0 || sum(abs(argu/l2n(argu))) <= sumabs) 
        return(0)
    lam1 <- 0
    lam2 <- max(argu) - 1e-05
    iter <- 1
    while (iter <= 15 && (lam2 - lam1) > (1e-04)) {
        su <- soft(argu, (lam1 + lam2)/2)
        if (sum(abs(su/l2n(su))) < sumabs) {
            lam2 <- (lam1 + lam2)/2
        } else {
            lam1 <- (lam1 + lam2)/2
        }
        iter <- iter + 1
    }
    return((lam1 + lam2)/2)
}
