#sets
process = 1:4
chemical = 1:6
scheme = 1:4
supplier = 1:4
plant = 1:3
customer = 1:4
intp = 1:20
intWW = zeros(length(intp))
for i in intp
  intWW[i] = 1.3^(i-1) - 1
end
J1 = 1:24
scenarios = 1:9
IJ = ((1,1),(1,6),(2,1),(2,6),(3,1),(3,2),(4,3),(4,4))
OJ = ((1,3),(2,3),(2,4),(3,3),(3,4),(4,5),(4,6))
PS = ((1,1),(2,1),(2,2),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2))
JM = ((1,1,3),(2,1,3),(2,2,4),(3,1,3),(3,2,4),(3,3,3),(3,4,4),(4,1,5),(4,2,5))
L = ((1,1,1),(1,1,6),(2,1,1),(2,2,1),(2,2,6),(3,1,1),(3,2,1),(3,3,2),(3,4,2),(4,1,6),(4,2,6))
Lbar = ((4,1,3),(4,2,4))
RJ = ((1,2),(1,4),(1,6),(2,1),(2,4),(2,6),(3,1),(3,2),(3,6),(4,1),(4,2),(4,4),(4,6))

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
basephi[1:4,  3] = 3.0
basephi[1:4, 5] = 300.0
QEU = zeros(3,4)
QEU[1:3, 1:3] = 1.50
QEU[1:3, 4] = 0.1
PUU = 100.0
FUU = 150.0

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

baseD = zeros(4,6)
baseD[1,3] = 100.0
baseD[1,5] = 5.0
baseD[2,3] = 50.0
baseD[2,5] = 3.0
baseD[3,3] = 150.0
baseD[3,5] = 2.5
baseD[4,3] = 80.0
baseD[4,5] = 2.0

baseS = zeros(4,6)
baseS[1,2] = 0.5
baseS[1,4] = 2.0
baseS[1,6] = 0.25
baseS[2,1] = 0.65
baseS[2,4] = 1.8
baseS[2,6] = 0.3
baseS[3,1] = 0.9
baseS[3,2] =1.0
baseS[3,6] = 0.25
baseS[4,1] = 0.75
baseS[4,2] = 0.75
baseS[4,4] = 1.9
baseS[4,6] = 0.3
betaS=[
      0.0     0.5      0.0         2.0          0.0         0.25
    0.65   0.0          0.0        1.8      0.0          0.30
    0.90        1.0       0.0       0.0       0.0          0.25
    0.75     0.75     0.0          1.9    0.0            0.30]
alphaRP = [
    20.0       10.0       15.0
    25.0       30.0       20.0
    30.0       15.0       20.0
    15.0       20.0       25.0]
betaRP =    [
      0.15         0.1       0.15
      0.12         0.2       0.05
      0.08        0.15       0.10
      0.13        0.08       0.20]
alphaPC =       [
         20.0          15.0        25.0        10.0
         15.0          25.0        30.0        35.0
         10.0          15.0        20.0        40.0]
betaPC =       [
       0.03         0.1      0.15      0.05
       0.20        0.15      0.05      0.20
       0.25        0.10      0.15      0.30]

#demand and price
baseprob = [0.25 0.5 0.25] 

#demand and price
#3 scenarios
scenarios=1:3
D = zeros(4,6,3)
D[1:length(customer), 1:length(chemical), 1] = baseD[1:length(customer), 1:length(chemical)] * 0.7
D[1:length(customer), 1:length(chemical), 2] = baseD[1:length(customer), 1:length(chemical)] 
D[1:length(customer), 1:length(chemical), 3] = baseD[1:length(customer), 1:length(chemical)] * 1.3
prob=baseprob
phi = basephi
#probability



