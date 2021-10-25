
# Simulations -------------------------------------------------------------

Data.generator.regression <- function(n,p,nb1,seedlist, NSIM, sdX=1){
  l1 <- list()
  for(q in 1:NSIM){
    set.seed(seedlist[q])  
    X          <- matrix( rnorm(n*p,sd=sdX),n,p)
    Xte        <- matrix( rnorm(n*p,sd=sdX),n,p)
    X          <- scale(X)
    Xte        <- scale(Xte)
    betat      <- c(rnorm(nb1,0,10),rnorm(nb1,0,1),rep(0,p-2*(nb1)))
    Y          <- X%*%betat + rnorm(n,0,1) + 5
    Yte        <- Xte%*%betat + rnorm(n,0,1) +5
    l1[[q]]         <- list(X=X,
                            Xte=Xte,
                            betat=betat,
                            Y=Y,
                            Yte=Yte)
  }
  return(l1)
}



Data.generator.mean <- function(n,nb1,seedlist,NSIM, sigma=1){
  l1 <- list()
  for(q in 1:NSIM){
    set.seed(seedlist[q])  
    betat      <- c(rnorm(nb1,0,10),rnorm(nb1,0,1),rep(0,n-2*(nb1)))
    Y          <- betat + rnorm(n,0,1)
    l1[[q]]         <- list(betat=betat,
                            Y=Y)
  }
  return(l1)
}
