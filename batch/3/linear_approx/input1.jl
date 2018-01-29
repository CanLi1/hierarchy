#sets
include("util.jl")
products=1:5
stages = 1:6
integer = 1:4
J1 = 1:(length(integer) * length(stages))
intp = 1:300


#parameters
alpha = ones(length(stages)) * 250.0
beta = ones(length(stages)) * 0.6
delta = 1330.0
lambda = ones(length(stages)) * 2000
VL = 300.0
VU = 2300.0
H = 5000

S =[
  7.9     2.0      5.2     4.9       6.1       4.2
   0.7     0.8      0.9     3.4       2.1       2.5
  0.7     2.6      1.6     3.6       3.2       2.9
  4.7     2.3      1.6     2.7       1.2       2.5
  1.2     3.6      2.4     4.5       1.6       2.1]

t = [
  6.4     4.7     8.3      3.9     2.1     1.2
  6.8     6.4     6.5      4.4     2.3     3.2
     1     6.3     5.4     11.9     5.7     6.2
   3.2       3     3.5      3.3     2.8     3.4
   2.1     2.5     4.2      3.6     3.7     2.2]  

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


#point of interpolation
int_tl_b = zeros(length(products), length(intp))
int_n_v = zeros(length(stages), length(intp))
int_ns = zeros(length(integer), length(stages), length(scenarios))
lb_tl_b = zeros(length(products))
ub_tl_b = zeros(length(products))

for i in products
	lb_tl_b[i] = max(t[i, 1:length(stages)]) / length(integer) /VU * max(S[i, 1:length(stages)])
	ub_tl_b[i] = max(t[i, 1:length(stages)])/VL * max(S[i, 1:length(stages)])
end

for i in products
	step = (ub_tl_b[i] - lb_tl_b[i]) / (length(intp) - 1) 
	for int in intp
		int_tl_b[i, int] = log(lb_tl_b[i] + (int-1) * step)
	end
end

lb_nv = VL^0.6
ub_nv = length(integer)*VU^0.6
step = (ub_nv-lb_nv) /(length(intp) - 1)
for j in stages
	for int in intp
		int_n_v[j, int] = log(lb_nv + (int - 1) * step)
	end
end



























