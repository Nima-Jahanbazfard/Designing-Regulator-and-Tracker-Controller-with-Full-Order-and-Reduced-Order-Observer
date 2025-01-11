%controller and observer


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

%two desierd poles for observer
p1=[-3+1i -3-1i -2 -3];
p2=[-4.5+1.5i -4.5-1.5i -3.5 -4.5];

%calculating L
AT=transpose(A);
CT=transpose(C);

%ackerman
Ka_1=acker(AT,CT,p1);
La1=1.*transpose(Ka_1);
Ka_2=acker(AT,CT,p2);
La2=1.*transpose(Ka_2);



%bass_gura
Kbg_1=Bass_Gura(AT,CT,p1);
Lbg1=transpose(Kbg_1);
Kbg_2=Bass_Gura(AT,CT,p2);
Lbg2=transpose(Kbg_2);



%equivalency
Ke_1=place(AT,CT,p1);
Le1=transpose(Ke_1);
Ke_2=place(AT,CT,p2);
Le2=transpose(Ke_2);


L1=Le1;
L2=Le2;



%Regulator



%determining the stability and cotrollability of the system
phi_c = ctrb(A,B);

if(rank(phi_c) == length(A))
display('phi_c is full rank so this system is controllable');
end

eigen = eig(A);
unstable_poles=eigen(real(eigen)>=0)
stable_poles=eigen(real(eigen)<0)
if isempty(unstable_poles)
    disp('this system is stable');
else 
    disp('this system is not stable');
end


%finding state feedback :
% 1-equivalency 2-Bass and Gura  3-Ackerman 4-canonical controller

%desired poles
desired_poles1 = [-1-0.3i -1+0.3i -0.7 -1];
desired_poles2 =  [-1.5-0.5i -1.5+0.5i -1.15 -1.5];

%equivalency
Ke1 = place(A, B, desired_poles1);
Ke2 = place(A, B, desired_poles2);

%Bass and Gura:the function of Bass_Gura is end of the code
Kbg1=Bass_Gura(A,B,desired_poles1);
Kbg2=Bass_Gura(A,B,desired_poles2);

%Ackerman
Ka1 = acker(A, B, desired_poles1);
Ka2 = acker(A, B, desired_poles2);

%canonical controller
Ac=[0 1 0 0;0 0 1 0;0 0 0 1; -1.062 -4 -3.5 0 ];
Bc=[0;0;0;1];
Cc=[-9.6 0.6 1 0];
Dc=0;
m=[0 0 0 1];
a1=m*(-1.*Ac);
s=tf('s');
delta_s_desired1=(s+1+0.3i)*(s+1-0.3i)*(s+0.7)*(s+1);
a2=[0.763 3.253 5.19 3.7];
delta_s_desired2=(s+1.5+0.5i)*(s+1.5-0.5i)*(s+1.15)*(s+1.5);
a3=[4.312  11.8 12.18 5.65];
Kc1=a2-a1;
Kc2=a3-a1;
phi_c_Ac=ctrb(Ac,Bc);
inv_phi_c=inv(phi_c);
Kcc1=Kc1*phi_c_Ac*inv_phi_c;
Kcc2=Kc2*phi_c_Ac*inv_phi_c;





%in total we have state feedbacks

%state feedback for desierd poles 1
K1=Ka1;
%state feedback for desierd poles 2
K2=Ka2;




%tracker static

%close loop for K1
Acl_1=A-B*K1;
%close loop for K2
Acl_2=A-B*K2;

%finding ua1 for first close loop system
ua1=inv((-1.*C)*inv(Acl_1)*B);

%finding ua2 for second close loop system
ua2=inv((-1.*C)*inv(Acl_2)*B);

%Tracking : Integral controller
%testing controllability of [B A;0 -C]
phi=[0 0 1 0 0;1 0 0 1 0;0.6 0 0 0 1;-13.1 -1.062 -4 -3.5 0;0 -1 0 0 0];
r=rank(phi);

if(r == 5)
display('[B A;0 -C] is full rank so this system is controllable');
display('so it is ok to designing trackig controller(integral)');
else
display('[B A;0 -C] is not full rank so this system is not controllable');

end
AI=[0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;-1.062 -4 -3.5 0 0;1 0 0 0 0];
BI=[4.520193743116709e-17; 1 ;0.6; -13.1; 0];
CI=[1 0 0 0 0];

dp1 = [-1-0.3i -1+0.3i -0.7 -1 -2];
dp2 =  [-1.5-0.5i -1.5+0.5i -1.15 -1.5 -1];
%calculating for both desired poles
KI_1=acker(AI,BI,dp1);
KI1_1=KI_1(:,1:4);
KI1_2=KI_1(:,5);
KI_2=acker(AI,BI,dp2);
KI2_1=KI_2(:,1:4);
KI2_2=KI_2(:,5);


%A_r for testing robustness
A_r=[0 0.8  0 0; 0 0 0.7 0 ;0 0 0 1;0 -3.5 -3 -1];
%B_r for testing robustness
B_r=[0;0.8;0.6;-11];
%B_r for testing robustness
C_r=[ 0.9  0.1 0 0];

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
