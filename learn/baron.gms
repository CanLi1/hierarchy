Variables
x
y;
Equations
e2;

e2 .. log(1+y) =g=x+10;
x.lo = 3; 
model a /all/;
option nlp = ipopt;
solve a using nlp minimizing y;
display e2.m;