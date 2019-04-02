function [ K,gamma,feas ] = Hinf_FilterK( SYS,GAMMA )
% HINF_OF return the gain for a luenberger Hinf filter 
% This type of oserver cannot account for polytopic uncertainties!
% D21==I

A=SYS.A{1}; B1=SYS.B1{1}; C1=SYS.C1{1}; C2=SYS.C2{1};

[nx,nw]= size(B1);
[ny,~]= size(C2);

% LMI toolbox solution
% setlmis([]);
% 
% if nargin==2
%     gam2=GAMMA^2;
%     J=lmivar(1,[1,1]);
%     lmiterm([1 2 2 0],-gam2);
%     lmiterm([1 3 3 0],-gam2);
%     MIN=0;
% else
%     gam2=lmivar(1,[1,1]);
%     lmiterm([1 2 2 gam2],-1,1);
%     lmiterm([1 3 3 gam2],-1,1);
%     MIN=1;
% end
% P=lmivar(1,[nx,1]); 
% Y=lmivar(2,[nx,ny]);
% % LMI 1:
% lmiterm([1 1 1 P],1,A,'s');
% lmiterm([1 1 1 Y],-1,C2,'s');
% lmiterm([1 1 1 0],C1'*C1);
% lmiterm([1 1 2 P],1,B1);
% lmiterm([1 1 3 Y],1,1);
% % LMI 2:
% lmiterm([2 1 1 P],-1,1);
% %lmiterm([2 1 1 0],-eps);
% % solve LMIs:
% lmisys=getlmis;
% options=[10^-4 1000 10^9 10 0];
% cobj=zeros(1+sum(1:nx)+ny*nx,1);
% cobj(1)=1;
% 
% if MIN==1
%     [~,xopt]=mincx(lmisys,cobj,options);
%     gamma2=dec2mat(lmisys,xopt,gam2);
%     gamma=sqrt(gamma2);
% else
%     [~,xopt]=feasp(lmisys,options);    
%     gamma=GAMMA;
% end
% 
% P=dec2mat(lmisys,xopt,P);
% Y=dec2mat(lmisys,xopt,Y);
% eig(P)
% K=P\Y;
% 
% evals = evallmi(lmisys,xopt);
% [lhs1,rhs1] = showlmi(evals,1) ;
% [lhs2,rhs2] = showlmi(evals,2);
% feas=all([eig(lhs1-rhs1) ; eig(lhs2-rhs2)]<0);

% yalmip solution:
if nargin==2
    gam2=GAMMA^2;
    MIN=0;
else
    gam2 = sdpvar(1);
    MIN=1;
end
P = sdpvar(nx,nx,'symmetric');
Y = sdpvar(nx,ny,'full');
lmi1 = [ P*A+A'*P-Y*C2-C2'*Y'+C1'*C1  P*B1           Y
         B1'*P                        -gam2*eye(nw)  zeros(nw,ny)
         Y'                           zeros(ny,nw)   -gam2*eye(ny) ];
     
constraints=[P>0,lmi1<0];
ops = sdpsettings ('solver','sdpt3');
if MIN==1
    solvesdp(constraints,gam2,ops);
else
    solvesdp(constraints,[],ops);
end
feas = all(checkset(constraints)>0);
gamma = sqrt(double(gam2));
K=double(P)\double(Y);

end

