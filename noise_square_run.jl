####Measure parameters
##k: the size of the support of measure to be considered
k = 5
##d: the degree of the moments to be considered
d = 5
#########################
#########################
##tol: the variance of the error distribution
tol = 1/100000

include("superres.jl")
include("noise_square.jl")
