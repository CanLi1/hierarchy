using JuMP
include("modelgen.jl")

a = now()
newmodel = spmodel()
b= now()
print(prate, "\n")
status = solve(newmodel)
print("status of the model is ", status, "\n")
print("Optimal objective value = ", getobjectivevalue(newmodel),"\n")


d = now()
newmodel2 = spmodel()
c = now()
status = solve(newmodel2)
print(typeof(status))
# print("status of the model is ", status, "\n")
# print("Optimal objective value = ", getobjectivevalue(newmodel2),"\n")

print(b-a, c-d)