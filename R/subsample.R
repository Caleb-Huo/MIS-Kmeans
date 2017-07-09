subsample <- function(aS) {
    wholeIndex = 1:nrow(aS)
    trainIndex = sample(wholeIndex, floor(length(wholeIndex)/2))
    testIndex = wholeIndex[-trainIndex]
    trainS = aS[trainIndex, ]
    testS = aS[testIndex, ]
    splitS = list(trainS = trainS, testS = testS)
    return(splitS)
}
