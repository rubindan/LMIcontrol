function [ K,gamma ] = HinfOFctrl( A,B1,B2,C1,C2,D12,D21,gamma_in )
%HINOFCTRL computes a stabilizing O.F H-inf controller
% A,B1,B2 are structs constructing the uncertain system
%   xdot = A{i}*x + B1{i}*w + B2{i}*u 
%   Z = C1*x + D12*u 
%   Y = C2*x + D21*w
% 
% Outputs the controller K(s) and attetuation level gamma
%
% HINOFCTRL requires preloading of YALMIP and SDPT3!

% Created: 26-May-2016 (Daniel Rubin)

% force a cell structure
if ~iscell(A), A={A}; end
if ~iscell(B1), B1={B1}; end
if ~iscell(B2), B2={B2}; end
%if ~iscell(C2), C2={C2}; end
%if ~iscell(D21), D21={D21}; end
N=max([length(A),length(B1),length(B2)]);
if ~isequal(length(A),length(B1),length(B2))
    if length(A)==1, [A{1:N}] = deal(A{1}); end
    if length(B1)==1, [B1{1:N}] = deal(B1{1}); end
    if length(B2)==1, [B2{1:N}] = deal(B2{1}); end
end

[nx,nu]=size(B2{1});
ny=size(C2,1);

X = sdpvar(nx); % n by n symmetric
Y = sdpvar(nx); % n by n symmetric
Ac = sdpvar(nx,nx,'full');
Bc = sdpvar(nx,ny,'full');
Cc= sdpvar(nu,nx,'full');
Dc= sdpvar(nu,ny,'full');

if nargin<8, 
    gam2 = sdpvar(1);
    feas=0;
else
    gam2 = gamma_in^2;
    feas=1; % solve the feasability problem!
end

constraints = [];
for k =1:N
    lmi1= nestedOFlmi( A{k},B1{k},B2{k},C1,C2,D12,D21,gam2,X,Y,Ac,Bc,Cc,Dc );
    constraints =[constraints lmi1<0];
end
lmi2=[X eye(nx) ; eye(nx) Y];
constraints= [constraints lmi2>0];
sdpopt = sdpsettings('solver','sdpt3');

if feas==0 %find minimum gamma
    DIAGNOSTIC = solvesdp(constraints,gam2,sdpopt);
else %solve the feasability problem
    DIAGNOSTIC = solvesdp(constraints,[],sdpopt);
end
disp(yalmiperror(DIAGNOSTIC.problem));
if ~all(checkset(constraints)>0)
    error('problem is infeasible!')
end
XC=double(X);
YC=double(Y);
gamma=sqrt(double(gam2));
M=eye(nx)-XC*YC;
DC= double(Dc);
CC= double(Cc)-DC*C2*YC;
BC= M \ ( double(Bc)-XC*B2{k}*DC );
AC= M \ ( double(Ac)-XC*A{k}*YC-XC*B2{k}*CC-M*BC*C2*YC-XC*B2{k}*DC*C2*YC );

K = ss(AC,BC,CC,DC);

end

function [ LMI ] = nestedOFlmi ( A,B1,B2,C1,C2,D12,D21,gam2,X,Y,Ac,Bc,Cc,Dc )
% OFLMI returns the first LMI for the Hinf O.F problem
[~,nw]=size(B1); [nz,~]=size(C1);

LMI = [A'*X+X*A+Bc*C2+C2'*Bc'  Ac+A'+C2'*Dc'*B2'        X*B1+Bc*D21    C1'+C2'*Dc'*D12'
       Ac'+A+B2*Dc*C2          Y'*A'+A*Y+B2*Cc+Cc'*B2'  B1+B2*Dc*D21   Y'*C1'+Cc'*D12'
       B1'*X+D21'*Bc'          B1'+D21'*Dc'*B2'        -gam2*eye(nw)   D21'*Dc'*D12'
       C1+D12*Dc*C2            C1*Y+D12*Cc              D12*Dc*D21     -eye(nz)         ];

end
