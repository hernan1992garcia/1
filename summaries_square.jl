#########################
##such that n+1 is the dimension of the points
n = 1
##########Domain
##########
roof = -1
ceil = 1

th_horiz = roof:0.01:ceil
th_vert = roof:0.01:ceil
###Domain set
@polyvar x[1:2]

dominio = [1-x[1]^2,1-x[2]^2]

mons_half = []
for i in 0:d
    append!(mons_half,monomials(x,i))
end
mons_half = convert(Array{typeof(mons_half[1])},mons_half)

mons = []
for i in 0:2*d
    append!(mons,monomials(x,i))
end
mons = convert(Array{typeof(mons[1])},mons)

mons_module = []
for i in 0:2*s
    append!(mons_module,monomials(x,i))
end
mons_module = convert(Array{typeof(mons_module[1])},mons_module)

############################################
lim_inf = [-1,-1]
lim_sup = [1,1]

mons_ort = mons
matrix_moments = zeros(length(mons_ort),length(mons_ort))

for i in 1:length(mons_ort)
    mono = mons_ort[i]
    matrix_moments[:,i] = [coefficient(mono,mons[k]) for k in 1:length(mons)]
end

mons_moments = Array{Any}(length(mons))

for i in 1:length(mons)
    mons_moments[i] = exponents(mons[i])
end

mons_moments_half = Array{Any}(length(mons_half))

for i in 1:length(mons_half)
    mons_moments_half[i] = exponents(mons_half[i])
end

###generate a measure
matrizCoef = matrix_moments
######################################

integ_vec = zeros(length(mons_ort))

for i in 1:length(mons_ort)
    integ_vec[i] = integ_poly(lim_inf,lim_sup,mons_ort[i])
end

delta = 0.0001
######################################

integ_v = clean_measurements_square(mons,mons_ort,delta,integ_vec,dominio)

(A,errr) = error_SOS_aprox(integ_vec,mons_ort,mons_half,matrizCoef)

#######################################
value_m = zeros(length(th_horiz),length(th_vert))

for u in 1:length(th_horiz)
    for v in 1:length(th_vert)
        punto = [th_horiz[u],th_vert[v]]
        ev_basis = evalua_polys_high([punto],mons_ort[1:length(mons_half)])
        value_m[u,v] = ((ev_basis)'*A*ev_basis)[1]
    end
end

pcolormesh(th_horiz,th_vert,log.(abs.((value_m)')))
colorbar()
