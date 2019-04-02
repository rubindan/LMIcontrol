
%% initialize problem 1:
clc
clear all
a=[ -0.5 0.5]; % a= alpha
A=@(n) [0 1+0.5*a(n) ; -1 -0.5];
B1 =[ -2 0 ; 1 0];
C1 =[1 0];
C2 =100*[ -1 1];
D21 =[0 0.8];
s=tf('s'); % define transfer function parameter

for k=1:2
   SYS.A{k}=A(k);
   SYS.B1{k}=B1;
   SYS.C1{k}=C1;
   SYS.C2{k}=C2;
   SYS.D21{k}=D21; 
end

[ FILTER,gamma,feas ] = Hinf_RobustFilter( SYS )