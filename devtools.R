WD <- '/Users/caleb/Box Sync/twowaySKM/MIS-Kmeans'
system(paste0("mkdir -p ", "'",WD, "'"))
setwd(WD)

if(F){
	install.packages("formatR")	
}

packageName <- "MISKmeans"
devtools::create(packageName) 

## licenses
setwd(WD)
setwd(packageName)

## put some code in the R package

## make the code neat
formatR::tidy_dir("R")

## prepare .Rd refer to the slides

## create .Rd R documentation
devtools::document() 

## useDynLib(MISKmeans) 
## // [[Rcpp::export]]


if(F){
  ## Exported data
  ## documenting data
  xxxx <- sample(1000)
  devtools::use_data(xxxx ,internal=TRUE)

  ## create test functions
  devtools::use_testthat()
  if(F){
    ## actual test function
    test_that("test if f function is correct", {
      expect_equal(f(1,1), 2)
    })
  }

  ## perform test
  devtools::test()
  
  
  ## browse all Vignettes
  browseVignettes()

  ## create our own Vignettes
  devtools::use_vignette("Algebra2")
  
}



## check the package
devtools::check()

## build the package
devtools::build()


## install the package
remove.packages("MISKmeans")
devtools::install()

install.packages("../MISKmeans_0.0.0.9000.tar.gz",repos=NULL,type="source")


## use the package
library(MISKmeans)

G <- 4
S <- 3
x <- NULL
for(s in 1:S){
  x[[s]] <- rnorm(G)
}

test(x,S,G)

if(F){
  browseVignettes("MISKmeans")
  browseVignettes()
}

help(package="MISKmeans")

if(F){
  nstart = 20
  ntrial = 1
  maxiter = 20
  lambda = 1/2
  method = "exhaustive"
  sampleSizeAdjust = FALSE
  wsPre = NULL
  silence = FALSE
  K<- 3
  wbounds <- 10
  awbound <- 10
  
  S = list(t(S1),t(S2), t(S2))
  res = MetaSpaKmeans(x=S,K=3,wbounds=10,lambda=2)
}

