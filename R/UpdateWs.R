UpdateWs <- function(x, Cs, l1bound, ratio, permatch) {
    a = ratio + permatch
    lam <- BinarySearch(a, l1bound)
    ws.unscaled <- soft(a, lam)
    return(ws.unscaled/l2n(ws.unscaled))
}
