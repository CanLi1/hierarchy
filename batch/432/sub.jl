using JuMP
using Mosek
using Ipopt
using KNITRO
using BARON
# function generate_sub(; djc=[[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15],[16],[17],[18],[19],[20],[21],[22],[23],[24]], Q=zeros(5), prob=0.0)
function generate_sub(; djc=[[1],[2],[3],[4],[5],[6],[7],[8],[9]], Q=zeros(5), prob=0.0)
    # s1 = Model(solver=KnitroSolver())
    # s1 = Model(solver=IpoptSolver())
    s1 = Model(solver=MosekSolver())
    #, KTR_PARAM_OUTLEV=1
	# s1 = Model(solver=BaronSolver())
	#KTR_PARAM_MAXTIMECPU=120.0
	ϵ=1e-5
    #sets for number of disjunctions
    disjunction=1:length(djc)

    #maximum number of disjuncts 
    max_djc = 1:512

	@NLparameter(s1, yfbar[int in  integer, j in stages]==0.0)
	@NLparameter(s1, nbar[j in stages]==0.0)
	@NLparameter(s1, vbar[j in stages]==0.0)

	s1[:yfbar] = yfbar
	s1[:nbar] = nbar
	s1[:vbar] = vbar

	@variable(s1, 0<=yf[int in integer, j in stages]<=1)
	@variable(s1, 0<=ys[j in J1]<=1)
	@variable(s1, n[j in stages]>=0)
	@variable(s1, log(VL)<=v[j in stages]<=log(VU))
	@variable(s1, ns[j in stages]>=0)
	@variable(s1, tl[i in products])
	@variable(s1, b[i in products])
	@variable(s1, L>=0)
	@variable(s1, firststagecost>=0)
	@variable(s1, secondstagecost>=0)

	@variable(s1, 0<=yff[int in integer, j in stages]<=1)
	@variable(s1, nn[j in stages])
	@variable(s1, log(VL)<=vv[j in stages]<=log(VU))

	#variable for disjunction 
	@variable(s1, 0<=dot_yf[int in integer, j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]<=1)
	@variable(s1, 0<=dot_y[j in J1, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]<=1)
	@variable(s1, dot_n[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
	@variable(s1, dot_v[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))])
	@variable(s1, dot_ns[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
	@variable(s1, dot_tl[i in products, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))])
	@variable(s1, dot_b[i in products, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))])
	@variable(s1, dot_L[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
	@variable(s1, dot_firststagecost[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
	@variable(s1, dot_secondstagecost[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
	@variable(s1, dot_λ[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
	
	#x=x̂
	@NLconstraint(s1, t01[int in integer, j in stages], yff[int, j]==yfbar[int, j])
	@NLconstraint(s1, t02[j in stages], nn[j]==nbar[j])
	@NLconstraint(s1, t03[j in stages], vv[j]==vbar[j])
	@constraint(s1, t1[int in integer, j in stages], yf[int, j] == yff[int, j])
	@constraint(s1, t2[j in stages], v[j] == vv[j])
	@constraint(s1, t3[j in stages], n[j] == nn[j])		

	# Ax+g(y) ⩽0, g2(y)⩽0
	@constraint(s1, c1[i in products, j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_v[j,d,k] >= log(S[i,j])*dot_λ[d,k] + dot_b[i,d,k])
	@constraint(s1, c2[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_ns[j,d,k] == sum(log(int) * dot_y[(int-1)*length(stages)+j,d,k] for int in integer))
	@constraint(s1, c3[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_ns[j,d,k] <= dot_n[j,d,k])
	@constraint(s1, c4[i in products, j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_ns[j,d,k] + dot_tl[i,d,k] >= log(t[i,j])*dot_λ[d,k])
	@NLconstraint(s1, c5[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], sum(Q[i] * (dot_λ[d,k]*(1-ϵ)+ϵ)*exp(dot_tl[i,d,k]/(dot_λ[d,k]*(1-ϵ)+ϵ)-dot_b[i,d,k]/(dot_λ[d,k]*(1-ϵ)+ϵ))- Q[i]*ϵ*(1-dot_λ[d,k]) for i in products) <= H*dot_λ[d,k] + dot_L[d,k])
	@constraint(s1, c6[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], sum(dot_y[(int-1)*length(stages)+j,d,k] for int in integer) == dot_λ[d,k])	
	@NLconstraint(s1, c7[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_firststagecost[d,k]>= prob * sum(alpha[j] * (dot_λ[d,k]*(1-ϵ)+ϵ)*exp(dot_n[j,d,k]/(dot_λ[d,k]*(1-ϵ)+ϵ) + beta[j] * dot_v[j,d,k]/(dot_λ[d,k]*(1-ϵ)+ϵ)) - alpha[j] * ϵ*(1-dot_λ[d,k]) for j in stages))
	@NLconstraint(s1, c8[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_secondstagecost[d,k]>=prob *( sum( lambda[j] *(dot_λ[d,k]*(1-ϵ)+ϵ) * exp(dot_ns[j,d,k]/(dot_λ[d,k]*(1-ϵ)+ϵ) ) -lambda[j]*ϵ*(1-dot_λ[d,k]) for j in stages)+ delta*dot_L[d,k]) )

    #x=∑ẋ, y=∑ẏ, ∑λ=1
    @constraint(s1, e1[int in integer, j in stages, d in disjunction], yf[int, j]==sum(dot_yf[int,j,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e2[j in J1, d in disjunction], ys[j] == sum(dot_y[j, d, k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e3[j in stages, d in disjunction], n[j] == sum(dot_n[j,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e4[j in stages, d in disjunction], v[j] == sum(dot_v[j,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e5[j in stages, d in disjunction], ns[j] == sum(dot_ns[j,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e6[i in products, d in disjunction], tl[i] == sum(dot_tl[i,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e7[i in products, d in disjunction], b[i] == sum(dot_b[i,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e8[d in disjunction], L == sum(dot_L[d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e9[d in disjunction], firststagecost == sum(dot_firststagecost[d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e10[d in disjunction], secondstagecost == sum(dot_secondstagecost[d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e11[d in disjunction], sum(dot_λ[d,k] for k in max_djc if k<= 2^(length(djc[d]))) == 1)

    #0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
    @constraint(s1, b1[int in integer, j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_yf[int,j,d,k]<= dot_λ[d,k])
    @constraint(s1, b2[j in J1, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_y[j, d,k] <= dot_λ[d,k])
    @constraint(s1, b3[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_n[j,d,k]<=dot_λ[d,k]*log(length(integer)))
    @constraint(s1, b4[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_v[j,d,k]>=dot_λ[d,k]*log(VL))
    @constraint(s1, b5[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_v[j,d,k]<=dot_λ[d,k]*log(VU))
    @constraint(s1, b6[j in stages, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_ns[j,d,k]<=dot_λ[d,k]*log(length(integer)))
    @constraint(s1, b7[i in products, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_tl[i,d,k]<=dot_λ[d,k]*log(max(t[i,:])))
    @constraint(s1, b8[i in products, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_tl[i,d,k]>=dot_λ[d,k]*log(max(t[i,:])/length(integer)))
    @constraint(s1, b9[i in products, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_b[i,d,k]<=dot_λ[d,k]*log(min(VU*inv_S[i,:])))
    @constraint(s1, b10[i in products, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_b[i,d,k]>=dot_λ[d,k]*log(min(VL*inv_S[i,:])))
    @constraint(s1, b11[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_L[d,k]<=dot_λ[d,k]*1.2e5)
    @constraint(s1, b12[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_firststagecost[d,k] <=prob * dot_λ[d,k] * 2.4e5)
    @constraint(s1, b13[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_secondstagecost[d,k] <= prob * dot_λ[d,k] * 2e7)

    #yⱼ=0 or 1
    index_d = 1
    for item in djc
        for index_k in 1:(2^length(item))
            seq_i = 1
            for index_i in item
                if (floor(Int, (index_k-1) / 2^(seq_i-1) ) % 2) == 1
                    @constraint(s1, dot_y[index_i, index_d, index_k] == dot_λ[index_d, index_k])
                else
                    @constraint(s1, dot_y[index_i, index_d, index_k] == 0)
                end
                seq_i = seq_i + 1
            end 
        end
        index_d = index_d + 1
    end

    @objective(s1, Min, firststagecost + secondstagecost)

    return s1 
end





