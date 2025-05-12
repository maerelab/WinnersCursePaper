#Fit spline and return prediction
splinPred = function(x, y, xOut, df){
    splin = smooth.spline(x, y, df = df)
    predict(splin, xOut)$y
}