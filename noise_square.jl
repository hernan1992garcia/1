
##n: such that K \subseteq R^{n+1}
n = 1
##########Domain
##########
roof = -1
ceil = 1

th_horiz = roof:0.05:ceil
th_vert = roof:0.05:ceil

#########################
######2s: degree of squared module
s = 6
#####Limits of integration
lim_inf = [roof,roof]
lim_sup = [ceil,ceil]
####Constructing a measure
indP_x = sample(1:length(th_horiz),k,replace = false)
indP_y = sample(1:length(th_vert),k,replace = false)

punt_supp = []
for i in 1:k
    append!(punt_supp,[[th_horiz[indP_x[i]],th_vert[indP_y[i]]]])
end
measureCoefs = (1/k)*ones(k)

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
mons_ort = ortonormalization(mons,lim_inf,lim_sup)
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
measureCoefs = (1/k)*ones(k)
matrizCoef = matrix_moments

#############Constructing b vector
integ_vec = evalua_polys_high(punt_supp,mons_ort)

Id = eye(length(measureCoefs))
for i in 1:length(measureCoefs)
    Id[i,i] = measureCoefs[i]
end

integ_vec = integ_vec*Id
b_vec = sum([integ_vec[:,i] for i in 1:k])

##noise
dist = Normal(0.0,tol)
noise = rand(dist,length(b_vec))

##b_noise: noised vector of measurements
b_noise = b_vec + noise
delta = norm(noise)

############################################
(integ_v,a_opti) = clean_measurements_square(mons,mons_ort,delta,b_noise,dominio)

(A,errr) = error_SOS_aprox(integ_v,mons_ort,mons_half,matrizCoef)
maximum(abs.(integ_v-b_noise))
#######################Plot
#######################
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


#for u in 1:length(punt_supp)
#    p = punt_supp[u]
#    plot(p[1],p[2],marker = "*",color=:red)
#end
