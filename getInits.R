getInits = function() { 

scale = 1.5
SDscale = 4

result = list()

result[["intercept"]] = sign(startingValues[["intercept" ]]) *
    runif(length(startingValues[["intercept" ]]),
       abs(startingValues[["intercept"]])/scale,
       scale * abs(startingValues[["intercept"]]))


result[["SDsite"]] = sqrt(runif(1,
       startingValues$vars[["site"]]/scale,
       startingValues$vars[["site"]]*scale))

result[["Rsite"]] = rnorm(length(startingValues[["Rsite"]]),
        startingValues[["Rsite"]], startingValues$vars[["site"]]/SDscale)

result[["phisite"]] = runif(1,
       startingValues$phi[["site"]]/scale,
       startingValues$phi[["site"]]*scale)


return(result)

}