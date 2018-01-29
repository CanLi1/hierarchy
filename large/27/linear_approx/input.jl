#sets
process = 1:4
chemical = 1:6
scheme = 1:4
supplier = 1:2
plant = 1:2
customer = 1:2
intp = 1:20
intWW = zeros(length(intp))
for i in intp
  intWW[i] = 1.3^(i-1) - 1
end
    #set of binary variable in the second stage
J1 = 1:8

IJ = ((1,1),(1,6),(2,1),(2,6),(3,1),(3,2),(4,3),(4,4))
OJ = ((1,3),(2,3),(2,4),(3,3),(3,4),(4,5),(4,6))
PS = ((1,1),(2,1),(2,2),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2))
JM = ((1,1,3),(2,1,3),(2,2,4),(3,1,3),(3,2,4),(3,3,3),(3,4,4),(4,1,5),(4,2,5))
L = ((1,1,1),(1,1,6),(2,1,1),(2,2,1),(2,2,6),(3,1,1),(3,2,1),(3,3,2),(3,4,2),(4,1,6),(4,2,6))
Lbar = ((4,1,3),(4,2,4))
RJ = ((1,2),(1,4),(1,6),(2,1),(2,4),(2,6))

#parameters
Q0=0.0
rho = 1.0
H = 1.0
betaC = [1.0 0.9 1.1 8]
alphaC = [20.0 10.0 15.0 30.0]
delta = zeros(4,4)
delta[1:3, 1:4] = 0.1
delta[4, 1:4] = 3.0
basephi = zeros(4,6)
basephi[1:4,  3] = 3.5
basephi[1:4, 5] = 100.0
QEU = zeros(3,4)
QEU[1:3, 1:3] = 0.75
QEU[1:3, 4] = 0.075
PUU = 75.0
FUU = 75.0

mu = zeros(4,4,6)
mu[1,1,1] = 1.05
mu[1,1,3] = -1.0
mu[1,1,6] = 0.03
mu[2,1,1] = 1.02
mu[2,1,3] = -1.0
mu[2,2,1] = 1.1
mu[2,2,4] = -1.0
mu[2,2,6] = 0.09
mu[3,1,1] = 1.1
mu[3,1,3] = -1.0
mu[3,2,1] = 1.2
mu[3,2,4] = -1.0
mu[3,3,2] = 1.08
mu[3,3,3] = -1.0
mu[3,4,2] = 1.05
mu[3,4,4] = -1.0
mu[4,1,3] = 1.2
mu[4,1,5] = -1.0
mu[4,1,6] = 1.0
mu[4,2,4] = 1.1
mu[4,2,5] = -1.0
mu[4,2,6] = 0.85

baseD = zeros(2,6)
baseD[1,3] = 100.0
baseD[1,5] = 5.0
baseD[2,3] = 50.0
baseD[2,5] = 3.0


betaS=[
      0.0     0.5      0.0         2.0          0.0         0.25
    0.65   0.0          0.0        1.8      0.0          0.30]

alphaRP = [
    20.0       10.0      
    25.0       30.0  ]     

betaRP =    [
      0.15         0.1      
      0.12         0.2      ]
alphaPC =       [
         20.0          15.0       
         15.0          25.0        ]
betaPC =       [
       0.08         0.1      
       0.20        0.15      
]

baseprob = [0.25 0.5 0.25] 

#demand and price
#3 scenarios
# scenarios=1:3
# D = zeros(2,6,3)
# D[1:2, 1:6, 1] = baseD[1:2,  1:6] * 0.7
# D[1:2, 1:6, 2] = baseD[1:2,  1:6] 
# D[1:2, 1:6, 3] = baseD[1:2,  1:6] * 1.3
# prob=baseprob
scenarios = 1:27       
prob = zeros(27)
D = zeros(2,6,27)
phi = zeros(2, 6,27)
ratios = [0.7 1.0 1.3]
for gen1 in 1:3
  for gen2 in 1:3
    for gen3 in 1:3
      temp_s = 3*(gen1-1) + gen2 + 9 * (gen3-1)
      D[1, 1:6, temp_s] = baseD[1, 1:6] *  ratios[gen1]
      D[2, 1:6, temp_s] = baseD[2, 1:6] *  ratios[gen2]
      prob[temp_s] = baseprob[gen1] * baseprob[gen2] * baseprob[gen3] 
      phi[1:2, 1:6, temp_s] = basephi[1:2, 1:6] * ratios[gen3]
    end
  end
end
# #81 scenarios
# scenarios = 1:81       
# prob = zeros(81)
# D = zeros(2,6,81)
# ratios = [0.7 1.0 1.3]
# for gen1 in 1:3
#   for gen2 in 1:3
#     for gen3 in 1:3
#       for gen4 in 1:3
#         temp_s = 3*(gen1-1) + gen2 + 9 *(gen3-1) + 27 * (gen4-1)
#         D[1, 3, temp_s] = baseD[1, 3] *  ratios[gen1]
#         D[1, 5, temp_s] = baseD[1, 5] * ratios[gen2]
#         D[2, 3, temp_s] = baseD[2, 3] *  ratios[gen3]
#         D[2, 5, temp_s] = baseD[2, 5] * ratios[gen4]
#         prob[temp_s] = baseprob[gen1] * baseprob[gen2] * baseprob[gen3] * baseprob[gen4]
#       end
#     end
#   end
# end

#probability



