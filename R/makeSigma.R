# Helper function to build covariance matrix
makeSigma = function(Sd, p){
    Sigma = matrix(Sd, p, p)
    diag(Sigma) = 1
    Sigma
}