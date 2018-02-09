sets
i /i1*i3/
a /a1*a2/
p /p1*p4/;
alias (i, j)
parameters
L(i) /i1 6, i2 6, i3 3/
H(i) /i1 6, i2 6, i3 3/
xbar(a) /a1 0, a2 0/
ybar(a) /a1 5, a2 30/
C(i, j)
base_Carea(a) /a1 800, a2 800/
M /100/
M2 /1500/
xlb(i) 
xub(i)
ylb(i)
yub(i);
C(i,j) = 0;
C('i1', 'i2') = 30;
C('i1', 'i3') = 24;
C('i2', 'i3') = 10;
xlb(i) = -10;
xub(i) = 10;
ylb(i) = -5;
yub(i) = 40;

sets
s /s1*s100/
sub /1*10/;

parameters
prob(s)
ratio(sub) /1 1.5, 2 1.4, 3 1.3, 4 1.2, 5 1.1, 6 1, 7 0.9, 8 0.8, 9 0.7, 10 0.6 /
Carea(a,s);
prob(s) = 1 / 100;
alias(sub, sub1);
loop(sub,
	loop(sub1,
		loop(s,
		if(ord(s) eq ((ord(sub)-1) * 10 + ord(sub1)),
	Carea('a1', s) = base_Carea('a1') * ratio(sub);
	Carea('a2', s) = base_Carea('a2') * ratio(sub1) ;
	);
	);
	);
	);





positive variables
delx(i,j), dely(i,j),rsqr(a,s);

variables
x(i), y(i), xs(i,s), ys(i,s), obj;

Binary variables
Z(i,j,p), W(i,a,s);

Equations
e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17;

e1(i,j)$(ord(i)<ord(j)) .. delx(i,j) =g= x(i) - x(j);
e2(i,j)$(ord(i)<ord(j)) .. delx(i,j) =g= x(j) - x(i);
e3(i,j)$(ord(i)<ord(j)) .. dely(i,j) =g= y(i) - y(j);
e4(i,j)$(ord(i)<ord(j)) .. dely(i,j) =g= y(j) - y(i);
e5(i,j)$(ord(i)<ord(j)) .. x(i) + L(i)/2 =l= x(j) - L(j) /2 + M*(1 - Z(i,j,'p1'));
e6(i,j)$(ord(i)<ord(j)) .. x(j) + L(j)/2 =l= x(i) - L(i) /2 + M*(1 - Z(i,j,'p2')) ;
e7(i,j)$(ord(i)<ord(j)) .. y(i) + H(i)/2 =l= y(j) - H(j) /2 + M*(1 - Z(i,j, 'p3'));
e8(i,j)$(ord(i)<ord(j)) .. y(j) + H(j)/2 =l= y(i) - H(i) /2 + M*(1 - Z(i,j,'p4'));
e9(i,j)$(ord(i)<ord(j)) .. sum(p, Z(i,j,p)) =e= 1;
e10(i,a,s) .. power((xs(i,s) - L(i)/2 - xbar(a)), 2) + power((ys(i,s) + H(i)/2 - ybar(a)), 2) =l= rsqr(a,s) + M2*(1-W(i,a,s));
e11(i,a,s) .. power((xs(i,s) - L(i)/2 - xbar(a)), 2) + power((ys(i,s) - H(i)/2 - ybar(a)), 2) =l= rsqr(a,s) + M2*(1-W(i,a,s));
e12(i,a,s) .. power((xs(i,s) + L(i)/2 - xbar(a)), 2) + power((ys(i,s) + H(i)/2 - ybar(a)), 2) =l= rsqr(a,s) + M2*(1-W(i,a,s));
e13(i,a,s) .. power((xs(i,s) + L(i)/2 - xbar(a)), 2) + power((ys(i,s) - H(i)/2 - ybar(a)), 2) =l= rsqr(a,s) + M2*(1-W(i,a,s));
e14(i,s) .. sum(a, W(i,a,s)) =e= 1;
e15(i,j,s)$(ord(i)<ord(j)) .. xs(i,s) - xs(j,s) =e= x(i) - x(j);
e16(i,j,s)$(ord(i)<ord(j)) .. ys(i,s) - ys(j,s) =e= y(i) - y(j);
e17 .. obj =e= sum((i,j)$(ord(i)<ord(j)), C(i,j) * (dely(i,j) + delx(i,j)) ) + sum(s, prob(s)* sum(a, rsqr(a,s) * Carea(a,s) ) );

model clay /all/;
clay.optfile=1;
options optca = 0;
options optcr =0.0001;
  OPTION LIMROW = 0;
option minlp=sbb;
clay.nodlim=1e8;  
OPTION LIMCOL = 0;
option reslim=1e4;
OPTION threads=12;
solve clay using minlp minimizing obj;





















