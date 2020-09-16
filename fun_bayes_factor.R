Bf <- function(se, # note that this is an error in the script and should be "se" 
                   #   (instead of "sd"); if we insert se instead of sd, 
                   # we obtain the exact same results as Dienes' online calculator 
                   # at www.lifesci.sussex.ac.uk/home/Zoltan_Dienes/inference/bayes_factor.swf 
               obtained, uniform, lower=0, upper=1, meanoftheory=0, sdtheory=1, tail=2)
{
  #test data can be found starting at plOO of Dienes (2008)
  #
  area <- 0 
  if(identical(uniform, 1)) {
    theta <- lower
    range <- upper - lower 
    incr <- range / 2000 
    for (A in -1000:1000) {
      theta <- theta + incr
      dist_theta <- 1 / range
      height <- dist_theta * dnorm(obtained, theta, se) 
      area <- area + height * incr
      } 
    } else{
      theta <- meanoftheory - 5 * sdtheory
      incr <- sdtheory / 200 
      for (A in -1000:1000){
        theta <- theta + incr
        dist_theta <- dnorm(theta, meanoftheory, sdtheory) 
        if(identical(tail, 1)){
          if (theta <= 0){ 
            dist_theta <- 0
          } else {
            dist_theta <- dist_theta * 2
            }
        }
        height <- dist_theta * dnorm(obtained, theta, se) 
        area <- area + height * incr
      }
    }
  LikelihoodTheory <- area
  Likelihoodnull <- dnorm(obtained, 0, se) 
  BayesFactor <- LikelihoodTheory / Likelihoodnull 
  ret <- list("LikelihoodTheory" = LikelihoodTheory, 
              "Likelihoodnull" = Likelihoodnull, 
              "BayesFactor" = BayesFactor)
  ret 
}

