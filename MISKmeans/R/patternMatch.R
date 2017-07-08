patternMatch <- function(x, Cs, ws, method = method, silence = FALSE) {
    if (method == "exhaustive") {
        matchExhaustive2(x, Cs, ws, silence = silence)
        return(matchExhaustive(x, Cs, ws, silence = silence))
    } else if (method == "linear") {
        return(matchLinear(x, Cs, ws, silence = silence))
    } else if (method == "MCMC") {
        return(matchMCMC(x, Cs, ws, silence = silence))
    } else {
        stop("Please specify method")
    }
}
