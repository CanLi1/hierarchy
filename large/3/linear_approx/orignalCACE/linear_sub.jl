using JuMP
using Gurobi
function generate_sub(; djc=[[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15],[16],[17],[18],[19],[20],[21],[22],[23],[24]], demand=zeros(4,6), price=zeros(4,6), xbar=zeros(3,4), Qbar=zeros(3,4), prob=0.0)
# function generate_sub(; djc=[[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12]], demand=zeros(4,6), price=zeros(4,6), xbar=zeros(3,4), Qbar=zeros(3,4), prob=0.0)
# function generate_sub(; djc=[[1],[2],[3],[4],[5],[6],[7],[8]], demand=zeros(2,6), price=zeros(2,6),  prob=0.0, Qbar=zeros(length(plant), length(process)), xbar=zeros(length(plant), length(process)))
	# s1 = Model(solver=MosekSolver(MSK_IPAR_NUM_THREADS=1, MSK_IPAR_INTPNT_MAX_ITERATIONS=1000000, MSK_DPAR_INTPNT_TOL_REL_GAP=1e-5, MSK_DPAR_INTPNT_NL_TOL_REL_GAP=1e-5))
	# s1 = Model(solver=IpoptSolver(tol=1e-6, constr_viol_tol=1e-3, compl_inf_tol=1e-3))
    # s1 = Model(solver=CplexSolver(CPX_PARAM_THREADS=1))
    s1 = Model(solver=GurobiSolver())
	ϵ=1e-5
    #sets for number of disjunctions
    disjunction=1:length(djc)

    #maximum number of disjuncts 
    max_djc = 1:256


 	#original variables. all the binary variables are expressed in single y 
    @variable(s1, PU[r in supplier, p in plant, j in chemical; (r,j) in RJ]>=0.0)
	@variable(s1, F[p in plant, c in customer, j in chemical]>=0.0)
	@variable(s1, QE[p in plant, i in process]>=0.0)
	@variable(s1, theta[p in plant, i in process, j in chemical, s in scheme; ((i,s) in PS && (i,s,j) in JM)]>=0.0)
	@variable(s1, WW[p in plant, i in process, j in chemical, s in scheme; ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)]>=0.0)
	@variable(s1, Slack[c in customer, j in chemical]>=0.0)
	@variable(s1, 0.0<=x[p in plant, i in process]<=1.0)
	@variable(s1, 0.0<=y[j in J1]<=1.0)

	#variable for disjunction
    @variable(s1, dot_PU[r in supplier, p in plant, j in chemical, d in disjunction, k in max_djc; (r,j) in RJ && k<= 2^(length(djc[d]))]>=0.0)
	@variable(s1, dot_F[p in plant, c in customer, j in chemical, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0.0)
	@variable(s1, dot_QE[p in plant, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0.0)
	@variable(s1, dot_theta[p in plant, i in process, j in chemical, s in scheme, d in disjunction, k in max_djc; k<= 2^(length(djc[d])) && ((i,s) in PS && (i,s,j) in JM)]>=0.0)
	@variable(s1, dot_WW[p in plant, i in process, j in chemical, s in scheme, d in disjunction, k in max_djc; k<= 2^(length(djc[d])) && ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar) ]>=0.0 )
	@variable(s1, dot_Slack[c in customer, j in chemical, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0.0)
	@variable(s1, 0.0<=dot_x[p in plant, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]<=1.0)
	@variable(s1, 0.0<=dot_y[j in J1, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]<=1.0)
	@variable(s1, dot_λ[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)


    # x=x̄
    @constraint(s1, t1[p in plant, i in process], QE[p,i]== Qbar[p,i])
    @constraint(s1, t2[p in plant, i in process], x[p,i]== xbar[p,i])

    # Ax+g(y) ⩽0, g2(y)⩽0
    @constraint(s1, c1[p in plant, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_QE[p,i,d,k] <= QEU[p,i] * dot_x[p,i,d,k])
	@constraint(s1, c3[p in plant, j in chemical, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], sum( dot_PU[r,p,j,d,k] for r in supplier if (r,j) in RJ) + sum(dot_WW[p,i,j,s,d,k] for i in process, s in scheme if ( (i,j) in OJ && (i,s) in PS && ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)) ) == sum(dot_F[p,c,j,d,k] for c in customer) + sum(dot_WW[p,i,j,s,d,k]  for i in process, s in scheme if ((i,j) in IJ && (i,s) in PS && ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar))))
	@constraint(s1, c4[p in plant, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], sum(dot_theta[p,i,j,s,d,k] for j in chemical, s in scheme if ((i,s,j) in JM && (i,s) in PS)) <= H * dot_QE[p,i,d,k] * 100.0 )
	@constraint(s1, c5[p in plant, i in process, j in chemical, s in scheme, d in disjunction, k in max_djc; (i,s) in PS && (i,s,j) in JM && k<= 2^(length(djc[d]))], dot_WW[p,i,j,s,d,k] == rho * dot_theta[p,i,j,s,d,k])	
	@constraint(s1, c7[p in plant, i in process, j in chemical, jj in chemical, s in scheme, d in disjunction, k in max_djc; (i,s,j) in L && (i,s) in PS && (i,s,jj) in JM && k<= 2^(length(djc[d]))], dot_WW[p,i,j,s,d,k] == mu[i,s,j] * dot_WW[p,i,jj,s,d,k])
	@constraint(s1, c8[p in plant, i in process, j in chemical, jj in chemical, s in scheme, int in intp, d in disjunction, k in max_djc; (i,s,j) in Lbar && (i,s) in PS && (i,s,jj) in JM && k<= 2^(length(djc[d]))], -dot_λ[d,k] * log(1+intWW[int]) -1/(1+intWW[int]) * (dot_WW[p,i,j,s,d,k]-intWW[int] * dot_λ[d,k]) + mu[i,s,j] * dot_WW[p,i,jj,s,d,k] <=0.0)
	@constraint(s1, c9[r in supplier, p in plant, j in chemical, d in disjunction, k in max_djc; (r,j) in RJ && k<= 2^(length(djc[d]))], dot_PU[r,p,j,d,k] <= PUU * dot_y[length(plant) * (r-1) + p,d,k])
	@constraint(s1, c10[p in plant, c in customer, j in chemical, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_F[p,c,j,d,k] <= FUU * dot_y[length(plant) * length(supplier) + length(customer) *(p-1) + c, d,k])
	@constraint(s1, c11[c in customer, j in chemical, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], sum(dot_F[p,c,j,d,k] for p in plant) + dot_Slack[c,j,d,k] == demand[c,j] * dot_λ[d,k])
	

    #x=∑ẋ, y=∑ẏ, ∑λ=1
    @constraint(s1, e1[r in supplier, p in plant, j in chemical, d in disjunction; (r,j) in RJ], PU[r,p,j] == sum(dot_PU[r,p,j,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e2[p in plant, c in customer, j in chemical, d in disjunction], F[p,c,j] == sum(dot_F[p,c,j,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e3[p in plant, i in process, d in disjunction], QE[p,i] == sum(dot_QE[p,i,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e4[p in plant, i in process, j in chemical, s in scheme, d in disjunction; ((i,s) in PS && (i,s,j) in JM)], theta[p,i,j,s] == sum(dot_theta[p,i,j,s,d,k] for k in max_djc if k<= 2^(length(djc[d])) ))
    @constraint(s1, e5[p in plant, i in process, j in chemical, s in scheme, d in disjunction; ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)], WW[p,i,j,s] == sum(dot_WW[p,i,j,s,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e6[c in customer, j in chemical, d in disjunction], Slack[c,j] == sum(dot_Slack[c,j,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e7[p in plant, i in process, d  in disjunction], x[p,i] == sum(dot_x[p,i,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e8[j in J1, d in disjunction], y[j] == sum(dot_y[j,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e9[d in disjunction], sum(dot_λ[d,k] for k in max_djc if k<= 2^(length(djc[d]))) == 1.0 )

    #0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
    @constraint(s1, b1[r in supplier, p in plant, j in chemical, d in disjunction, k in max_djc; (r,j) in RJ && k<= 2^(length(djc[d]))], dot_PU[r,p,j,d,k] <= PUU * dot_λ[d,k])
    @constraint(s1, b2[p in plant, c in customer, j in chemical, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_F[p,c,j,d,k] <= FUU * dot_λ[d,k])
    @constraint(s1, b3[p in plant, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_QE[p,i,d,k] <= QEU[p,i] * dot_λ[d,k])
    @constraint(s1, b4[p in plant, i in process, j in chemical, s in scheme, d in disjunction, k in max_djc; k<= 2^(length(djc[d])) && ((i,s) in PS && (i,s,j) in JM)], dot_theta[p,i,j,s,d,k] <= 100.0 * QEU[p,i] * dot_λ[d,k])
    @constraint(s1, b5[p in plant, i in process, j in chemical, s in scheme, d in disjunction, k in max_djc; k<= 2^(length(djc[d])) && ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)], dot_WW[p,i,j,s,d,k] <= 1.2*75.0 * dot_λ[d,k] )
    @constraint(s1, b6[p in plant, s in scheme, d in disjunction, k in max_djc; k<= 2^(length(djc[d])) && ((4,s,5) in JM || (4,s,5) in L || (4,s,5) in Lbar)], dot_WW[p,4,5,s,d,k] <= 10.0* dot_λ[d,k] )
    @constraint(s1, b7[c in customer, j in chemical, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_Slack[c,j,d,k] <= demand[c,j]*dot_λ[d,k])
    @constraint(s1, b8[p in plant, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_x[p,i,d,k]<= dot_λ[d,k])
    @constraint(s1, b9[j in J1, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_y[j,d,k] <= dot_λ[d,k])

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

    #objective
    @objective(s1, Min, prob * (sum(betaC[i] * QE[p,i] * 100.0 + alphaC[i] * x[p,i] for p in plant, i in process) + sum(delta[i,s] * rho * theta[p,i,j,s] for p in plant, i in process, s in scheme, j in chemical if ((i,s) in PS && (i,s,j) in JM)) + sum((betaS[r,j] + betaRP[r,p])  * PU[r,p,j] for p in plant, j in chemical, r in supplier if (r,j) in RJ) + sum(alphaRP[r,p] * y[length(plant) * (r-1) + p] for r in supplier, p in plant) + sum(alphaPC[p,c] * y[length(plant) * length(supplier) + length(customer) *(p-1) + c] for p in plant, c in customer) + sum(betaPC[p,c] * F[p,c,j] for p in plant, c in customer, j in chemical) + sum(price[c,j] * Slack[c,j] for c in customer, j in chemical) ) )
    
	return s1 
end