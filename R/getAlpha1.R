getAlpha1 = function(mat){
    corMat = cor(mat)
    mean(corMat[upper.tri(corMat)])
}