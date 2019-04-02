function [ K,gam ] = dHinfSFGain( A,B1,B2,C1,D11,D12,gamma )
%DHINFSFGAIN computes a DISCRETE stabilizing SF H-inf controller
% A,B1,B2 are structs constructing the uncertain system
%   x(k+1) = A{i}*x(k) + B2{i}*u(k) + B1{i}*w(k)
% C1,D11,D12 are matrices determining the objective
%   z(k) = C1*x(k) + D12*u(k) + D11*w(k)
% 
% DHINFSFGAIN requires preloading YALMIP and SDPT3!
%
% Created: 2-October-2015 (Daniel Rubin)
% Last update: 21-August-2016 (Daniel Rubin)

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

[nx,nu]=size(B2{1});

Q = sdpvar(nx);
Y = sdpvar(nu,nx);

if nargin<7, 
    gam2 = sdpvar(1);
    feas=0;
else
    gam2 = gamma^2;
    feas=1; % solve the feasability problem!
end
constraints = Q>0;
for k =1:length(A)
    lmi1= nestedDSFlmi(A{k},B1{k},B2{k},C1,D11,D12,gam2,Q,Y );
    constraints =[constraints lmi1<=0];
end
ops = sdpsettings('solver','sedumi');

if feas==0 %find minimum gamma
    DIAGNOSTIC = solvesdp(constraints,gam2,ops);
else %solve the feasability problem
    DIAGNOSTIC = solvesdp(constraints,[],ops);
end
if DIAGNOSTIC.problem~=0
    error(yalmiperror(DIAGNOSTIC.problem));
end
% if ~all(checkset(constraints)>0)
%     error('problem is infeasible!')
% end
Q=double(Q);
Y=double(Y);
gam=sqrt(double(gam2));
K=Y/Q;

% check quadratic stability: (added: 21-August-2016)
P = sdpvar(nx);
constraints = [P>0];
for k=1:length(A)
    Ac = A{k}+B2{k}*K;
    LMI_k= [ Ac'*P*Ac-P ];
    constraints =[ constraints LMI_k<0 ];
end
ops = sdpsettings ('solver','sedumi');
DIAGNOSTIC = solvesdp(constraints,[],ops);
if all(checkset(constraints)>0)
    fprintf('\n ** Closed loop is quadraticaly stable! ** \n\n')
else
    fprintf('\n ** Closed loop is NOT quadraticaly stable! ** \n\n')
end

end

function [ LMI ] = nestedDSFlmi ( A,B1,B2,C1,D11,D12,gam2,Q,Y )
% SFLMI returns the first LMI for the Hinf S.F problem
[nz,nx]=size(C1);
[~,nw]=size(B1);
%[~,nu]=size(B2);

% In Eq. (88) from lecture notes repalce A,C1 with A+B2*K, C1+D12*K. Than
% take Y=K*Q and obtain: 
LMI=[ -Q           Q*A'+Y'*B2'  zeros(nx,nw)  Q*C1'+Y'*D12'
      A*Q+B2*Y     -Q           B1            zeros(nx,nz)
      zeros(nw,nx) B1'          -gam2*eye(nw) D11'
      C1*Q+D12*Y   zeros(nz,nx) D11           -eye(nz)      ]; 

% LMI = [ -Q           A*Q+B2*Y    zeros(nx,nz)   B1
%         Q*A'+Y'*B2'  -Q           Q*C1'+Y'*D12' zeros(nx,nw)
%         zeros(nz,nx) C1*Q+D12*Y   -gam2*eye(nz) D11 
%         B1'          zeros(nw,nx) D11'          -eye(nw)    ];

end
