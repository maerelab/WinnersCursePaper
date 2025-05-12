makeFolds = function(nfolds, data, cvType){
    folds = sample(nfolds, nrow(data), replace = TRUE)
    while(any(table(folds)<=1)){
        folds = makeFolds(nfolds, data, cvType)
    }
    return(folds)
}
