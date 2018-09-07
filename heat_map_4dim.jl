
##n: such that K \subseteq R^{n+1}
n = 3
##########Domain
##########
roof = 0
ceil = 1

th_1 = roof:0.1:ceil
th_2 = roof:0.1:ceil
th_3 = roof:0.1:ceil
th_4 = roof:0.1:ceil

th_mesh = []
for i1 in th_1
    for i2 in th_2
        for i3 in th_3
            for i4 in th_4
                append!(th_mesh,[[i1, i2,i3,i4]])
            end
        end
    end
end

M = zeros(size_k,size_d)
############################################
@polyvar x[1:n+1]
###################
d = 1
mons = []
for i in 0:d
    append!(mons,monomials(x,i))
end
mons = convert(Array{typeof(mons[1])},mons)

###################

for d in 1:size_d

    mons_half = []
    for i in 0:convert(Int64,floor(d/2))
        append!(mons_half,monomials(x,i))
    end

    mons = []
    for i in 0:d
        append!(mons,monomials(x,i))
    end
    mons = convert(Array{typeof(mons[1])},mons)

    matrizCoef = eye(length(mons))

    for k in 1:size_k

        measureCoefs = (1/k)*ones(k)

        prom_perc_recov = zeros(rep)

        for cont in 1:rep
        indP_x = sample(1:length(th_1),k)
        indP_y = sample(1:length(th_2),k)
        indP_z = sample(1:length(th_3),k)
        indP_w = sample(1:length(th_4),k)

        punt_supp = []
        for i in 1:k
            append!(punt_supp,[[th_1[indP_x[i]],th_2[indP_y[i]],th_3[indP_z[i]],th_4[indP_w[i]]]])
        end

          #############Constructing b vector
          integ_vec = evalua_polys(punt_supp,mons)

          Id = eye(length(measureCoefs))
          for i in 1:length(measureCoefs)
              Id[i,i] = measureCoefs[i]
          end

          integ_vec = integ_vec*Id
          b_vec = sum([integ_vec[:,i] for i in 1:k])
tic()
          (A,errr) = error_SOS_aprox(b_vec,mons,mons_half,matrizCoef)
         ########################
          ev_basis = evalua_polys(collect(th_mesh),mons[1:length(mons_half)])
          ev_basis_sup = evalua_polys(punt_supp,mons[1:length(mons_half)])

          zetas_sup = [log((ev_basis_sup[:,i])'*A*ev_basis_sup[:,i]+1) for i in 1:length(punt_supp)]

          is_or_not = 1*(zetas_sup.<= 0.0000001)
######################################################
        perc_recover = sum(is_or_not)/k
toc()
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
