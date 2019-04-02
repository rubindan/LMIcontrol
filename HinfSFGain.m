function [ K,gam ] = HinfSFGain( A,B1,B2,C1,D12,gamma )
%HINFSFGAIN computes a stabilizing S.F H-inf controller
% A,B1,B2,C1,D12 are struct constructing the uncertain system
%   xdot = A{i}*x + B2{i}*u + B1{i}*w
%   Z    = C*x + D12*u
%
% gamma:    if gamma is empty --> return K which minimizes gamma (default)
%           if gamma is a positive real scalar --> return corresponding K
%           if gamma equals -1 --> solve the feasibility problem

if nargin<6, gamma=[]; end

% force a cell structure
if ~iscell(A), A={A}; end
if ~iscell(B1), B1={B1}; end
if ~iscell(B2), B2={B2}; end
if ~isequal(length(A),length(B1),length(B2))
    N=max([length(A),length(B1),length(B2)]);
    if length(A)==1, [A{1:N}] = deal(A{1}); end
    if length(B1)==1, [B1{1:N}] = deal(B1{1}); end
    if length(B2)==1, [B2{1:N}] = deal(B2{1}); end
end

R=D12'*D12;
feas=0; 

[nx,nw]=size(B1{1});
nu=size(B2{1},2);
nz=size(C1,1);
X = sdpvar(nx);
Y = sdpvar(nu,nx);

if isempty(gamma), 
    gam2 = sdpvar(1);
else
    if gamma==-1,
        gam2 = sdpvar(1);
        feas=1;
    else
        gam2 = gamma^2;
    end
end

constraints=[X >0];
for k =1:length(A)
    lmi1{k}= nestedSFlmi(A{k},B1{k},B2{k},C1,R,gam2,X,Y );
    constraints =[ constraints lmi1{k}<0];
end
ops = sdpsettings ('solver','sdpt3','verbose',2,'sdpt3.maxit',500);
%ops = sdpsettings ('solver','sedumi','verbose',2);
if ~feas 
    DIAGNOSTIC = solvesdp(constraints,gam2,ops); %find minimum gamma
else 
    DIAGNOSTIC = solvesdp(constraints,[],ops); %solve the feasability problem
end
disp(yalmiperror(DIAGNOSTIC.problem));
if ~all(checkset(constraints)>-10^-4)
    error('problem detected!')
end

gam= sqrt ( double ( gam2 ));
X= double (X);
Y= double (Y);
K=Y/X;

end

function [ LMI ] = nestedSFlmi ( A,B1,B2,C1,R,gam2,X,Y )
% SFLMI returns the first LMI for the Hinf S.F problem
[nz,nx]=size(C1);
[~,nw]=size(B1);
[~,nu]=size(B2);

LMI= [ A*X+X*A'+B2*Y+Y'*B2' B1            Y'           X*C1'
       B1'                  -gam2*eye(nw) zeros(nw,nu) zeros(nw,nz)
       Y                    zeros(nu,nw)  -inv(R)      zeros(nu,nz)
       C1*X                 zeros(nz,nw)  zeros(nz,nu) -eye(nz)     ];
end
