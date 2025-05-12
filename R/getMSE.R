#Extract MSE from simulation object
getMSE = function(x){
    sapply(x, function(y) y$MSE)
}
