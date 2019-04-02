function [ FILTER,gamma,feas ] = Hinf_RobustFilter( SYS,GAMMA )
% HINF_ROBUSTFILTER return the general filter stracture for a robust Hinf
% filter.  
%   SYS - structure with system matrices A,B1,C1,C2,D21
%   GAMMA - if not defined, solver will try to minimize gamma 

if nargin==2
    gam2=GAMMA^2;
    MIN=0;
else
    gam2 = sdpvar(1);
    MIN=1;
end
N=length(SYS.A);
[nx,nw]= size(SYS.B1{1});
[nz,~]= size(SYS.C1{1});
[nu,~]= size(SYS.C2{1});
X = sdpvar(nx);
R = sdpvar(nx);
Af = sdpvar(nx,nx,'full');
Bf = sdpvar(nx,nu,'full');
Cf = sdpvar(nz,nx,'full');
Dc = sdpvar(nz,nu,'full') ;

lmi2 = [X R ; R R];
constraints=lmi2>0;
for k=1:N
    A=SYS.A{k}; B1=SYS.B1{k}; C1=SYS.C1{k}; C2=SYS.C2{k}; D21=SYS.D21{k};
    lmi1 = Hinf_RobustFilterLMI1( A,B1,C1,C2,D21,gam2,X,R,Af,Bf,Cf,Dc );
    constraints=[lmi1<0,constraints];
end
ops = sdpsettings ('solver','sdpt3');
if MIN==1
    solvesdp(constraints,gam2,ops);
else
    solvesdp(constraints,[],ops);
end
feas = all(checkset(constraints)>0);
FILTER.X = double(X);
FILTER.R = double(R);
FILTER.Af = double(Af);
FILTER.Bf = double(Bf);
FILTER.Cf = double(Cf);
FILTER.Dc = double(Dc);
gamma = sqrt(double(gam2));
% Gf=Cf*inv(s*(R-X)-Af)*Bf+Dc;
% Gf= minreal (Gf);
end

