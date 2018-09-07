
###points of set K
th = -1:0.001:1
########limits of integration
lim_inf = [-1]
lim_sup = [1]
####Constructing  measure
indP = sample(1:length(th),k,replace = false)
supp_pun = th[indP]
measureCoefs = (1/k)*ones(k)
######2s: degree of squared module
s = 11
##n: dimension of the sphere
n = 0
###Domain set
@polyvar x
dominio = [1-x^2]
####
mons_half = []
for i in 0:d
    append!(mons_half,monomials(x,i))
end

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
########################
########################

##mons_ort: orthonormal basis
mons_ort = ortonormalization(mons,lim_inf,lim_sup)

matrix_moments = zeros(length(mons_ort),length(mons_ort))
matrizCoef = matrix_moments
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

#############Constructing b vector
integ_vec = evalua_polys(supp_pun,mons_ort)

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
#################################
#################################
(integ_v,a_opti) = clean_measurements_line(mons_module,mons_ort,delta,b_noise,dominio)

(A,errr) = error_SOS_aprox(integ_v,mons_ort,mons_half,matrizCoef)

ev_basis = evalua_polys(collect(th),mons_ort[1:length(mons_half)])
ev_basis_sup = evalua_polys(supp_pun,mons_ort[1:length(mons_half)])

zetas = [log(abs((ev_basis[:,i])'*A*ev_basis[:,i])) for i in 1:length(th)]
zetas_sup = [log(abs((ev_basis_sup[:,i])'*A*ev_basis_sup[:,i])) for i in 1:length(supp_pun)]

####################Plot
#####################
plot(th,zetas,color =:blue)

for i in 1:length(supp_pun)
    y_line = minimum(zetas):0.01:maximum(zetas)
    plot(supp_pun[i]*ones(length(y_line)),y_line,linestyle = ":", color = :black)
    plot(supp_pun[i],minimum(zetas),marker = "*",color = :red)
end
