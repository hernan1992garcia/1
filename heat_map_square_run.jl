####Plot parameters
##size_d: maximum degree considered. Initial degree will be d = 1
size_d = 17
##size_k: maximum size support considered. Initial size support will be k = 1
size_k = 20
##rep: number of samples for each pair (k,d)
rep = 100

include("superres.jl")
include("heat_map_square.jl")
