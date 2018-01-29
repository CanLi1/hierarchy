#sets
products=1:5
stages = 1:3
integer = 1:3
J1 = 1:(length(integer) * length(stages))
#parameters
alpha = ones(length(stages)) * 250.0
beta = ones(length(stages)) * 0.6
delta = 230.0
lambda = ones(length(stages)) * 2000
VL = 300.0
VU = 2300.0
H = 5000

S =[
  7.9     2.0      5.2     
   0.7     0.8      0.9     
  0.7     2.6      1.6     
  4.7     2.3      1.6     
  1.2     3.6      2.4    ]

t = [
  6.4     4.7     8.3      
  6.8     6.4     6.5      
     1     6.3     5.4    
   3.2       3     3.5      
   2.1     2.5     4.2   ]  

inv_S = zeros(length(products), length(stages))
for i in products
  for j in stages
    inv_S[i,j] = 1.0/S[i,j]
  end
end
baseQ = zeros(length(products))
baseQ[1] = 250000 
baseQ[2]= 150000
baseQ[3]= 180000
baseQ[4]= 160000
baseQ[5]= 120000

baseprob = [0.3 0.4 0.3]
ratio =[1.1 1 0.9]

#3 scenarios
scenarios = 1:3
prob = baseprob
Q = zeros(length(products), 3)
Q[1:length(products), 1] = baseQ[1:length(products)] * ratio[1]
Q[1:length(products), 2] = baseQ[1:length(products)] * ratio[2]
Q[1:length(products), 3] = baseQ[1:length(products)] * ratio[3]

