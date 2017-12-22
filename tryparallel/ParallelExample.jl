using JuMP, CPLEX

type Data
   D::Array{Float64, 2}
   C::Array{Float64, 1}
   T::Array{Float64, 2}
   S::Array{Float64, 1}
   P::Array{Float64, 1}
end

function generateinstance(totalSuppliers::Int64, totalClients::Int64, totalScenarios::Int64)
   #generate a random instance of the problem
   D = rand(totalClients,totalScenarios)
   C = 100*rand(totalSuppliers)
   T = 10*rand(totalSuppliers, totalClients)
   S = ones(totalSuppliers)*100.0
   P = ones(totalScenarios).*1/totalScenarios
   return Data(D, C, T, S, P)
end

function generatemodel(ins::Data)
   #generate JuMP model
   Suppliers = 1:size(ins.C,1) #captures the total of suppliers
   Clients = 1:size(ins.D,1) #captures the total of clients
   Scenarios = 1:size(ins.P,1) #captures the total of scenarios

   m = Model(solver = CplexSolver(OutputFlag = 0))

   @variables(m, begin
      y[Suppliers], Bin
      x[Suppliers, Clients, Scenarios] >= 0
   end)

   @constraints(m, begin
      demand[j in Clients, s in Scenarios], sum(x[i,j,s] for i in Suppliers) >= ins.D[j,s]
      supply[i in Suppliers, s in Scenarios], sum(x[i,j,s] for j in Clients) <= ins.S[i]*y[i]
   end)

   @objective(m, Min,
      sum(ins.C[i]*y[i] for i in Suppliers) +
      sum(ins.P[s]*(sum(ins.T[i,j]*x[i,j,s] for i in Suppliers, j in Clients)) for s in Scenarios));

   return m
end

function generatemodels(totalSuppliers::Int64, totalClients::Int64, totalScenarios::Int64, totalModels::Int64)
   # generate an array of models of same size
   modelArray = Array(JuMP.Model,totalModels)

   for i = 1:totalModels
      modelArray[i] =
         generatemodel(generateinstance(totalSuppliers::Int64, totalClients::Int64, totalScenarios::Int64))
   end

   return modelArray
end

function solvemodel(m::JuMP.Model)
   # solve and print obj. function value for each problem
   solve(m)
   println(getobjectivevalue(m))
end

function generate_then_solvemodel(totalSuppliers::Int64, totalClients::Int64, totalScenarios::Int64)
   m = generatemodel(generateinstance(totalSuppliers, totalClients, totalScenarios))

   # solve and print obj. function value for each problem
   solve(m)
   println(getobjectivevalue(m))

   return "solved by worker $(myid())"
end