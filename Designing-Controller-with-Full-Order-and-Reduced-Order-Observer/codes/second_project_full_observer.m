
%determining the stability and cotrollability of the system
clc;
clear all;

num=[0 0 1 0.6 -9.6];
den=[1 0 3.5 4 1.062];

A=[0 1 0 0;0 0 1 0;0 0 0 1;-1.062 -4 -3.5 0];
B=inv([1 0 0 0;0 1 0 0;3.5 0 1 0;4 3.5 0 1])*[0;1;0.6;-9.6];
C=[1 0 0 0];

phi_o=obsv(A,C);
n = rank(phi_o);
if (n==size(A))
     disp('this system is observabel');
else
     disp('this system is not observabel');
end

%two desierd poles
p1=[-3+1i -3-1i -2 -3];
p2=[-4.5+1.5i -4.5-1.5i -3.5 -4.5];

%calculating L
AT=transpose(A);
CT=transpose(C);

%ackerman
Ka1=acker(AT,CT,p1);
La1=1.*transpose(Ka1);

Ka2=acker(AT,CT,p2);
La2=1.*transpose(Ka2);


%bass_gura
Kbg1=Bass_Gura(AT,CT,p1);
Lbg1=transpose(Kbg1);

Kbg2=Bass_Gura(AT,CT,p2);
Lbg2=transpose(Kbg2);


%equivalency
Ke1=place(AT,CT,p1);
Le1=transpose(Ke1);

Ke2=place(AT,CT,p2);
Le2=transpose(Ke2);
%in total we have
L1=Le1;
L2=Le2;




function k = Bass_Gura(A,B,pd)
phi_c = ctrb(A,B);
alpha = poly(pd);
alpha = alpha(1,2:end);
n = length(A);
e = eig(A);
a = poly(e);
a = a(1,2:end);
si=eye(n);
for i = 2:n
si = si+diag(a(i-1)*ones(1,n-i+1),i-1);
end
k = (alpha -a)*inv(si)*inv(phi_c);
end







