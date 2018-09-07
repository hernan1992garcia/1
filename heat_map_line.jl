
##n+1 = dimension
n = 0
## ceil, roof: such that K = [roof,ceil]
roof = 0.0
ceil = 1.0
th = roof:0.005:ceil

#################
M = zeros(size_k,size_d)
######################

###Domain set
@polyvar x

###beggining
d = 1
mons = []
for i in 0:d
    append!(mons,monomials(x,i))
end
mons = convert(Array{DynamicPolynomials.Monomial{true}},mons)
########

for d in 1:size_d
    for k in 1:size_k

        prom_perc_recov = zeros(rep)
        for cont in 1:rep

          indP = sample(1:length(th),k,replace = false)
          punt_supp = th[indP]

          mons_half = []
          for i in 0:convert(Int64,floor(d/2))
              append!(mons_half,monomials(x,i))
          end

          mons = []
          for i in 0:d
              append!(mons,monomials(x,i))
          end
          mons = convert(Array{DynamicPolynomials.Monomial{true}},mons)

          measureCoefs = (1/k)*ones(k)
          matrizCoef = eye(length(mons))

          #############Constructing b vector
          integ_vec = evalua_polys(punt_supp,mons)

          Id = eye(length(measureCoefs))
          for i in 1:length(measureCoefs)
              Id[i,i] = measureCoefs[i]
          end

          integ_vec = integ_vec*Id
          b_vec = sum([integ_vec[:,i] for i in 1:k])

          (A,errr) = error_SOS_aprox(b_vec,mons,mons_half,matrizCoef)
          ########################
          ev_basis = evalua_polys(collect(th),mons[1:length(mons_half)])
          ev_basis_sup = evalua_polys(punt_supp,mons[1:length(mons_half)])

          zetas_sup = [log((ev_basis_sup[:,i])'*A*ev_basis_sup[:,i]+1) for i in 1:length(punt_supp)]
          is_or_not = 1*(zetas_sup.<= 0.0000001)
######################################################
        perc_recover = sum(is_or_not)/k

        prom_perc_recov[cont] = perc_recover
        print(string("========*****ITERACION","(k,d)=",string(k),",",string(d)),"sample:",string(cont),"perc:",string(perc_recover))
    end

   close()
    M[k,d] = mean(prom_perc_recov)

   pcolormesh(M)
    colorbar()
   xlabel("d")
   ylabel("k")

    end
end
