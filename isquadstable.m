function [ STABLE,P ] = isquadstable( A,isdiscrete )
%ISQUADSTABLE check system's quadratic stability
% 
% ISQUADSTABLE(A) returns 1 if A is quadratically stable, and 0 if not. A
% ia a cell array containning all uncertain cases of the systems matrix A
%
% ISQUADSTABLE(A,ISDISCRETE) add an option to select betwin continous-time 
% (ISDISCRETE=0, default) and discrete-time (ISDISCRETE=1). 
%      
% ISQUADSTABLE requires preloading of YALMIP and SDPT3!

% Created: 23-May-2017 (Daniel Rubin)

if nargin<2, isdiscrete=[]; end
if isempty(isdiscrete), isdiscrete=0; end

% force a cell structure
if ~iscell(A), A={A}; end

[nx,~]=size(A{1});

% check quadratic stability: 
P = sdpvar(nx);
constraints = [P>0];
for k=1:length(A)
    if isdiscrete
        LMI_k= [ A{k}'*P*A{k}-P ]; % discrete time
    else
        LMI_k= [ P*A{k}+A{k}'*P ]; % continous time
    end
    constraints =[ constraints LMI_k < 0 ];
end
ops = sdpsettings('solver','sdpt3','verbose',2);
DIAGNOSTIC=solvesdp(constraints,[],ops);
disp(yalmiperror(DIAGNOSTIC.problem));

if all(checkset(constraints)>0)
    %fprintf('\n ** Closed loop is quadraticaly stable! ** \n\n')
    STABLE=1;
else
    STABLE=0;
end

P=double(P);

end