clc;
clear all;

num=[0 0 1 0.6 -9.6];
den=[1 0 3.5 4 1.062];

%make the system to observability realization
Ao=[0 1 0 0;0 0 1 0; 0 0 0 1;-1.062 -4 -3.5 0];
Bo=inv([1 0 0 0;0 1 0 0;3.5 0 1 0;4 3.5 0 1])*[0;1;0.6;-9.6];
Co=[1 0 0 0];
p=[1 0 0 0;0 0 0 1;0 1 0 0;0 0 1 0];
pinv=inv(p);
At=p*Ao*pinv;
Bt=p*Bo;
Ct=Co*pinv;
A11=At(1,1);
A12=At(1,2:4);
A21=At(2:4,1);
A22=At(2:4,2:4);
B1=Bt(1,1);
B2=Bt(2:4,1);

phi_o_reduced=obsv(A22,A12);
if (size(Ao)-1==rank(phi_o_reduced))
    disp('reduced system is observable');
else
    disp('reduced system is not onbservable');
end
%two desierd poles for reduced observer

pr1=[-3+1i -3-1i -2];
pr2=[-4.5+1.5i -4.5-1.5i -3.5];

k_reduced1=acker(transpose(A22),transpose(A12),pr1);
L_reduced1=transpose(k_reduced1);

k_reduced2=acker(transpose(A22),transpose(A12),pr2);
L_reduced2=transpose(k_reduced2);

%designing for first poles
A_1=A21-(L_reduced1*A11);
A_2=A22-(L_reduced1*A12);
B_1=B2-(L_reduced1*B1);
Q=inv(p);
Q1=Q(:,1);
Q2=Q(:,2:4);

%designing for second poles
A_1_1=A21-(L_reduced2*A11);
A_2_2=A22-(L_reduced2*A12);
B_1_1=B2-(L_reduced2*B1);
Q=inv(p);


%regulator
%desired poles
desired_poles1 = [-1-0.3i -1+0.3i -0.7 -1];
desired_poles2 =  [-1.5-0.5i -1.5+0.5i -1.15 -1.5];

%equivalency
Ker1 = place(Ao, Bo, desired_poles1)
Ker2 = place(Ao, Bo, desired_poles2)

%Bass and Gura:the function of Bass_Gura is end of the code
Kbgr1=Bass_Gura(Ao,Bo,desired_poles1)
Kbgr2=Bass_Gura(Ao,Bo,desired_poles2)

%Ackerman
Kar1 = acker(Ao, Bo, desired_poles1)
Kar2 = acker(Ao, Bo, desired_poles2)
%K1 & K2
K1=Kar1;
K2=Kar2;

%tracker static

%close loop for K1
Acl_1=Ao-Bo*K1;
%close loop for K2
Acl_2=Ao-Bo*K2;

%finding ua1 for first close loop system
ua1=inv((-1.*Co)*inv(Acl_1)*Bo)

%finding ua2 for second close loop system
ua2=inv((-1.*Co)*inv(Acl_2)*Bo)

%Tracking : Integral controller
AI=[Ao(1,:) 0; Ao(2,:) 0;Ao(3,:) 0;Ao(4,:) 0 ;Co 0];
BI=[Bo;0];
CI=[Co, 0];

%desierd poles
dp1 = [-1-0.3i -1+0.3i -0.7 -1 -2];
dp2 =  [-1.5-0.5i -1.5+0.5i -1.15 -1.5 -1];

%calculating for both desired poles
KI_1=acker(AI,BI,dp1)
KI1_1=KI_1(:,1:4);
KI1_2=KI_1(:,5);
KI_2=acker(AI,BI,dp2)
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
