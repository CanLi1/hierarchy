sets
i /i1*i3/
a /a1*a2/
p /p1*p4/
k /k1*k4/;
alias (i, j, ii, jj), (a,aa);
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
s /s1*s3/;

parameters
prob(s)/s1 0.3, s2 0.4, s3 0.3/
ratio(s) /s1 1.3, s2 1, s3 0.7/
Carea(a,s);

alias(s, s1);
loop(s1,
	Carea('a1', s1) = base_Carea('a1') * ratio(s1);
	Carea('a2', s1) = base_Carea('a2');
	);



positive variables
delx(i,j), dely(i,j),rsqr(a,s), dot_rsqr(a,s,i,aa);

variables
x(i), y(i), xs(i,s), ys(i,s), obj, dot_x(i, ii, jj, p), dot_y(i,ii,jj,p)
dot_xs(i,s,a), dot_ys(i,s,a);

Binary variables
Z(i,j,p), W(i,a,s);

equations
e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20,e21,e22,e23,e24,e25,e26,e27,e28,e29,e30,e31;

 e1(i,j)$(ord(i)<ord(j)) .. dot_x(i, i, j, 'p1') + L(i)/2 * Z(i,j, 'p1')  =l= dot_x(j, i, j, 'p1') - L(j) /2 * Z(i,j,'p1');
 e2(i,j)$(ord(i)<ord(j)) ..  dot_x(j, i, j, 'p2') + L(j)/2 * Z(i,j,'p2')  =l= dot_x(i, i,j,'p2') - L(i) /2 *Z(i,j,'p2');
 e3(i,j)$(ord(i)<ord(j)) ..  dot_y(i,i,j,'p3') + H(i)/2 *Z(i,j,'p3')  =l= dot_y(j,i,j,'p3') - H(j) /2 *Z(i,j,'p3');
e4(i,j)$(ord(i)<ord(j)) ..  dot_y(j,i,j,'p4') + H(j)/2*Z(i,j,'p4')  =l= dot_y(i,i,j,'p4') - H(i) /2 *Z(i,j,'p4');
 e5(i,j)$(ord(i)<ord(j)) ..  sum(p, Z(i,j,p)) =e= 1 ;

e6(i,ii,jj)$(ord(ii)<ord(jj)) .. x(i) =e= sum(p, dot_x(i, ii, jj, p) ) ;
e7(i,ii,jj)$(ord(ii)<ord(jj)) ..	 y(i) =e= sum(p, dot_y(i, ii, jj, p) );
e8(i,ii,jj,p)$(ord(ii)<ord(jj)) ..	 dot_x(i, ii, jj, p) =l= Z(ii,jj,p) * 10;
e9(i,ii,jj,p)$(ord(ii)<ord(jj)) ..	  dot_x(i, ii, jj, p)=g= Z(ii,jj,p) * (-10);
e10(i,ii,jj,p)$(ord(ii)<ord(jj)) ..	  dot_y(i, ii, jj, p) =l= Z(ii,jj,p) * 40;
e11(i,ii,jj,p)$(ord(ii)<ord(jj)) ..	  dot_y(i, ii, jj, p)=g= Z(ii,jj,p) * (-5);


e12(i,j)$(ord(i)<ord(j)) .. delx(i,j) =g= x(i) - x(j);
e13(i,j)$(ord(i)<ord(j)) ..   delx(i,j) =g= x(j) - x(i);
e14(i,j)$(ord(i)<ord(j)) ..   dely(i,j) =g= y(i) - y(j);
e15(i,j)$(ord(i)<ord(j)) ..   dely(i,j) =g= y(j) - y(i);

e16(i,j,s)$(ord(i)<ord(j)) ..   xs(i,s) - xs(j,s) =e= x(i) - x(j);
e17(i,j,s)$(ord(i)<ord(j)) ..    ys(i,s) - ys(j,s) =e= y(i) - y(j);



e18(i,a,s) .. 	 ((1-1e-6)*W(i,a,s) + 1e-6)*power((dot_xs(i,s,a)/((1-1e-6)*W(i,a,s) + 1e-6) - L(i)/2 - xbar(a)),2) - 1e-6  * power((-L(i)/2-xbar(a)),2)*(1-W(i,a,s)) + ((1-1e-6)*W(i,a,s) + 1e-6) * power((dot_ys(i,s,a)/((1-1e-6)*W(i,a,s) + 1e-6) + H(i)/2 - ybar(a)),2) - 1e-6 *power((H(i)/2 - ybar(a)),2)*(1-W(i,a,s))  =l= dot_rsqr(a,s,i,a) ;
e19(i,a,s) ..  	 ((1-1e-6)*W(i,a,s) + 1e-6)*power((dot_xs(i,s,a)/((1-1e-6)*W(i,a,s) + 1e-6) - L(i)/2 - xbar(a)),2) - 1e-6  * power((-L(i)/2-xbar(a)),2)*(1-W(i,a,s)) + ((1-1e-6)*W(i,a,s) + 1e-6) * power((dot_ys(i,s,a)/((1-1e-6)*W(i,a,s) + 1e-6) - H(i)/2 - ybar(a)),2) - 1e-6 *power((-H(i)/2 - ybar(a)),2)*(1-W(i,a,s))  =l= dot_rsqr(a,s,i,a) ;
e20(i,a,s) ..  	 ((1-1e-6)*W(i,a,s) + 1e-6)*power((dot_xs(i,s,a)/((1-1e-6)*W(i,a,s) + 1e-6) + L(i)/2 - xbar(a)),2) - 1e-6  * power((+L(i)/2-xbar(a)),2)*(1-W(i,a,s)) + ((1-1e-6)*W(i,a,s) + 1e-6) * power((dot_ys(i,s,a)/((1-1e-6)*W(i,a,s) + 1e-6) + H(i)/2 - ybar(a)),2) - 1e-6 *power((H(i)/2 - ybar(a)),2)*(1-W(i,a,s))  =l= dot_rsqr(a,s,i,a) ;
e21(i,a,s) ..  	 ((1-1e-6)*W(i,a,s) + 1e-6)*power((dot_xs(i,s,a)/((1-1e-6)*W(i,a,s) + 1e-6) + L(i)/2 - xbar(a)),2) - 1e-6  * power((+L(i)/2-xbar(a)),2)*(1-W(i,a,s)) + ((1-1e-6)*W(i,a,s) + 1e-6) * power((dot_ys(i,s,a)/((1-1e-6)*W(i,a,s) + 1e-6) - H(i)/2 - ybar(a)),2) - 1e-6 *power((-H(i)/2 - ybar(a)),2)*(1-W(i,a,s))  =l= dot_rsqr(a,s,i,a) ;
e22(i,s) ..  sum(a, dot_xs(i,s,a) ) =e= xs(i,s);
e23(i,s) ..  sum(a, dot_ys(i,s,a) ) =e= ys(i,s);
e24(a,s,i) ..	 rsqr(a,s) =e= sum(aa, dot_rsqr(a,s,i,aa) );
e25(i,s) .. 	 sum(a, W(i,a,s) ) =e= 1;
e26(i,s,a) .. 	 dot_xs(i,s,a)  =l= 10*W(i,a,s);
e27(i,s,a) .. 	 dot_xs(i,s,a) =g= (-10)*W(i,a,s);
e28(i,s,a) .. 	 dot_ys(i,s,a)  =l= 40*W(i,a,s);
e29(i,s,a) .. 	 dot_ys(i,s,a) =g= (-5)*W(i,a,s);
e30(a,s,i,aa) ..	 dot_rsqr(a,s,i,aa)  =l= 500 * W(i,aa,s);


e31 .. 	obj =e= sum((i,j)$(ord(i)<ord(j)), C(i,j) * (dely(i,j) + delx(i,j)) ) + sum(s, prob(s) * sum(a, rsqr(a,s) * Carea(a,s) ) );

model clay /all/;
clay.optfile=1;
options optca = 0;
options optcr =0.0001;
clay.nodlim=1e8;
option minlp=sbb;
  OPTION LIMROW = 0;
OPTION LIMCOL = 0;
option reslim=1e4;
OPTION threads=12;
solve clay using minlp minimizing obj; 

























