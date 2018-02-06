using JuMP
using Ipopt
using Mosek
using KNITRO 
using Pajarito
using CPLEX
#sets
rectangles = 1:3
areas = 1:2
pos = 1:4
#parameters
L = [6 6 3]
H = [6 6 3]
xbar = [0 0]
ybar = [5 30]
C = zeros(3,3)
C[1,2] = 30
C[1,3] =24
C[2,3] = 10
base_Carea = [800 800]
#big M
M = 100
M2 = 1500
xlb = ones(length(rectangles)) * (-10)
xub = ones(length(rectangles)) * 10
ylb = ones(length(rectangles)) * (-5)
yub = ones(length(rectangles)) * 40
#scenarios 
scenarios = 1:3
prob = [0.3 0.4 0.3]
Carea = zeros(length(areas), length(scenarios))
ratio = [1.3 1 0.7]
# ratio = [1 1 1]
for s in scenarios
	Carea[1, s] = base_Carea[1] * ratio[s]
	Carea[2, s] = base_Carea[2]
end