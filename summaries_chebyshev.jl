th = -1:0.001:1
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
lim_inf = [-1]
lim_sup = [1]

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

matrizCoef = matrix_moments

#############Constructing b vector-\frac{1}{\sqrt{1-x^2}}
##########################
integ_vec_mas = hquadrature(length(mons_ort), (w,v) -> v[:] = ev_mons_mas_cheb(w),0,1)[1]
integ_vec_menos = hquadrature(length(mons_ort), (w,v) -> v[:] = ev_mons_menos_cheb(w),0,1)[1]

integ_vec = integ_vec_mas + integ_vec_menos
delta = 0.0001
#########################################

(integ_v,a_opti) = clean_measurements_line(mons_module,mons_ort,delta,integ_vec,dominio)

(A,errr) = error_SOS_aprox(integ_v,mons_ort,mons_half,matrizCoef)

ev_basis = evalua_polys(collect(th),mons_ort[1:length(mons_half)])

zetas = [log(abs((ev_basis[:,i])'*A*ev_basis[:,i])) for i in 1:length(th)]
#####################################
plot(th,zetas,color =:blue)

cheb = chebyshev_zeros(d)

for i in 1:length(cheb)
    y_line = minimum(zetas):0.01:maximum(zetas)
    plot(cheb[i]*ones(length(y_line)),y_line,linestyle = ":", color = :black)
    plot(cheb[i],minimum(zetas),marker = "o",color = :red)
end
