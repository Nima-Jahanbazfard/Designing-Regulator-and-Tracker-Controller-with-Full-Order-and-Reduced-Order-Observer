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