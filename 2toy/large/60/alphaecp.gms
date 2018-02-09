Sets
s /1*60/;
alias (s, s1)
parameter
prob(s),d(s),unit_length;
unit_length = 3 / 60;
prob(s) = 1 / 60;
d('1') = 0;
loop(s1$(ord(s1)>1),
	d(s1) = d(s1-1) + unit_length;
	);

Positive variables
p1,p2,q1(s),q2(s),delta(s);

Binary variables
x1,x2,y1(s),y2(s);

variables
obj; 
Equations
e1,e2,e3,e4,e5,e6,e7,e8,e9,e10;

e1 .. p1 =l= 4 * x1;
e2 .. p2 =l= 2 * x2;
e3(s) .. power((q1(s) -3),2) + power((q2(s) - 2),2) =l= 1 + 16*(1-y1(s))  ;
e4(s) .. power((q1(s)-1),2) + power((q2(s)),2) =l=1+16*y1(s);
e5(s) .. power((q1(s)),2) + power((q2(s) - 1),2) =l= 1 + 16*(1-y2(s));
e6(s) .. power((q1(s)-4),2) + power((q2(s)-1),2) =l=1+16*y2(s);
e7(s) .. q1(s)=l= p1;
e8(s) .. q2(s)=l= p2;
e9(s) .. q1(s) + q2(s) + delta(s) =g= d(s);
e10 .. obj =e= sum(s, prob(s)* (q1(s)-12*q2(s) + 3*y1(s) -3*y2(s) + 100 * delta(s)) ) + 3 * x1 +3 * x2 + p1 + p2;
model toy /all/;
toy.optfile=1;
option optcr = 0.0001;
option optca = 0;
option reslim = 1e4;
option threads = 12;
option minlp = alphaecp;
solve toy using minlp minimizing obj;
