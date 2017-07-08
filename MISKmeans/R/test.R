##' @export
test <- function(x, S, G) {
    x <- as.double(x)
    S <- as.integer(S)
    G <- as.integer(G)
    .C("test_R", x, S, G)
}
