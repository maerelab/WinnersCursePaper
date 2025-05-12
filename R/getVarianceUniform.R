#Find the variance from the mixture of uniforms that define the prior for ash
getVarianceUniform = function(ashObj){
    with(ashObj$fitted, sum(pi*(b-a)^2/12))
}