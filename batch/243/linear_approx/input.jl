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
delta = 230
lambda = ones(length(stages)) * 2000
VL = 1300.0
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
scenarios = 1:243
prob = zeros(length(scenarios))
Q = zeros(length(products), 243)
for s1 in 1:3
  for s2 in 1:3
    for s3 in 1:3
      for s4 in 1:3
        for s5 in 1:3
          temp_s = 9 * (s1-1) + 3 * (s2-1) + s3 + 27 *(s4-1) + 81 * (s5-1)
          Q[1, temp_s] = baseQ[1] * ratio[s1]
          Q[2, temp_s] = baseQ[2] * ratio[s2]
          Q[3, temp_s] = baseQ[3] * ratio[s3]
          Q[4, temp_s] = baseQ[4] * ratio[s4]
          Q[5, temp_s] = baseQ[5] * ratio[s5]
          prob[temp_s] = baseprob[s1] * baseprob[s2] * baseprob[s3] * baseprob[s4] * baseprob[s5]
        end
      end
    end
  end
end


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



























