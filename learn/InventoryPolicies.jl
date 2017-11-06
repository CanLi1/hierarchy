module InventoryPolicies

using JuMP
using Gurobi

export model

include("InvPolicyData.jl")

for (k,t) in keys(Demand)
    Demand[k,t] == 0 && (Demand[k,t] = 10)
end

function model(policy, RC)
    m = Model(solver=GurobiSolver(MIPGap=0.002, Threads=4,OutputFlag=0))

@variable(m, supply[i in plants, j in warehouses, t in periods] >=0, upperbound=MaxSupply[i])
@variable(m, flow[j in warehouses, k in customers, t in periods] >= 0, upperbound=Demand[k,t])
@variable(m, inv[j in warehouses, t in (periods[1]-1):periods[end]] >=0, upperbound=MaxInv[j])
#@variable(m, z[i in plants, j in warehouses,t in periods], Bin)

for j in warehouses
    setupperbound(inv[j,periods[1]-1], InitialInv[j])
    setlowerbound(inv[j,periods[1]-1], InitialInv[j])
    setupperbound(inv[j,periods[end]], InitialInv[j])
    setlowerbound(inv[j,periods[end]], InitialInv[j])
end

    for (i,j,t) in keys(supply)
        if t <= L
            setupperbound(supply[i,j,t],0)
        end
    end

@constraint(m, demand[k in customers, t in periods], sum(flow[j,k,t] for j in warehouses) == Demand[k,t])
@constraint(m, invbal[j in warehouses, t in periods], inv[j,t] == inv[j,t-1] + sum(supply[i,j,t] for i in plants) - sum(flow[j,k,t] for k in customers))
#@constraint(m, semicontL[i in plants, j in warehouses, t in periods], z[i,j,t]*MinSupply[i] <= supply[i,j,t])
#@constraint(m, semicontR[i in plants, j in warehouses, t in periods], supply[i,j,t] <= z[i,j,t]*MaxSupply[i])

@expression(m, transportCostIJ, sum(TransportCost[i,j]*supply[i,j,t] for (i,j,t) in keys(supply) ))

@expression(m, transportCostIJ, sum(TransportCost[i,j]*supply[i,j,t] for (i,j,t) in keys(supply) ))
@expression(m, transportCostJK, sum(TransportCost[j,k]*flow[j,k,t] for (j,k,t) in keys(flow)))
#@expression(m, replenishCost, sum(ReplenishCost[i]*z[i,j,t] for (i,j,t) in keys(z) ))
@expression(m, inventoryCost, sum(HoldingCost*inv[j,t] for (j,t) in keys(inv)  )  + sum(HoldingCost*supply[i,j,t] for (i,j,t) in keys(supply)))
#@expression(m, inventoryCost, sum(HoldingCost*supply[i,j,t] for (i,j,t) in keys(supply)))
if RC
    @objective(m, Min, transportCostIJ + transportCostJK + replenishCost + 0.5*inventoryCost  )
else
    @objective(m, Min, transportCostIJ + transportCostJK + 0.05*inventoryCost)
end



if policy == :rQ
    @variable(m, y[j in warehouses, t in periods], Bin)
    @variable(m, r[j in warehouses] >= 0, upperbound=MaxInv[j])
    @variable(m, Q[j in warehouses] >= 0, upperbound=MaxInv[j])

    b = 0
    @constraint(m, qL1[j in warehouses, t in periods], sum(supply[i,j,t] for i in plants) >= Q[j] -b - MaxInv[j]*(1 - y[j,t]) )
    @constraint(m, qL2[j in warehouses, t in periods], sum(supply[i,j,t] for i in plants) <= Q[j] +b + MaxInv[j]*(1 - y[j,t]) )
    @constraint(m, qL[j in warehouses, t in periods], sum(supply[i,j,t] for i in plants) <= MaxInv[j]*y[j,t] )
    @constraint(m, minr1[j in warehouses, t in periods; t-L>=periods[1]], inv[j,t-L] <= r[j] +b + MaxInv[j]*(1 - y[j,t]) )
    @constraint(m, minr2[j in warehouses, t in periods; t-L>=periods[1]], inv[j,t-L] >= r[j] -b - MaxInv[j]*(1 - y[j,t]) )
end

if policy == :sS
    @variable(m, s[j in warehouses] >= 0, upperbound=MaxInv[j])
    @variable(m, S[j in warehouses] >= 0, upperbound=MaxInv[j])
    @variable(m, y[j in warehouses, t in periods], Bin)
    @variable(m, rp[j in warehouses, n in patterns], Bin)
    @variable(m, ra[j in warehouses, t in periods], Bin)
    @variable(m, rr[j in warehouses, t in periods], Bin)

    @constraint(m, allowedpolicy[j in warehouses, t in periods], sum(rp[j,n] for n in patterns if allowmatrix[n,t]>0) == ra[j,t])
    @constraint(m, onepattern[j in warehouses], sum(rp[j,n] for n in patterns) == 1)
    @constraint(m, allowedrepenishment[j in warehouses, t in periods], sum(supply[i,j,t]  for i  in plants) <= MaxInv[j]*y[j,t] )
    @constraint(m, inventoryconditionA[j in warehouses, t in periods; t-L>=periods[1]], sum(supply[i,j,t] for i in plants)  <= S[j] - inv[j,t-L] + MaxInv[j]*(1-y[j,t]) )
    @constraint(m, inventoryconditionB[j in warehouses, t in periods; t-L>=periods[1]], sum(supply[i,j,t] for i in plants)  >= S[j] - inv[j,t-L] - MaxInv[j]*(1-y[j,t]) )
    @constraint(m, inventoryconditionB[j in warehouses, t in periods; t-L>=periods[1]], sum(supply[i,j,t] for i in plants)  >= S[j] - inv[j,t-L] - MaxInv[j]*(1-y[j,t]) )
    @constraint(m, replenishrequiredA[j in warehouses, t in periods; t-L>=periods[1]], -MaxInv[j]*rr[j,t] + 0.001 <= inv[j,t-L] - s[j] )
    @constraint(m, replenishrequiredB[j in warehouses, t in periods; t-L>=periods[1]],  inv[j,t-L] - s[j] <= MaxInv[j]*(1- rr[j,t]) )
    @constraint(m, logicA[j in warehouses, t in periods], ra[j,t] + rr[j,t] <= y[j,t] + 1 )
    @constraint(m, logicB[j in warehouses, t in periods], y[j,t] <= ra[j,t] )
    @constraint(m, logicC[j in warehouses, t in periods], y[j,t] <= rr[j,t] )
    m.obj += sum(ra[j,t] for (j,t) in keys(ra))

end

return m
end


end