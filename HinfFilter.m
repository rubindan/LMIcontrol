function [ LMI ] = HinfFilter (  A,B1,C1,C2,D21,gam2,X,R,Af,Bf,Cf,Dc  )
%OFLMI Summary of this function goes here
%   Detailed explanation goes here
[~,nw]=size(B1); 
[nz,~]=size(C1);
LMI = [A'*X+X*A+Bf*C2+C2'*Bf'  X*A+Af+Bf*C2+A'*R  X*B1+Bf*D21   C1'-C2'*Dc'
       A'*X+Af'+C2'*Bf'+R*A    A'*R+R*A           R*B1          C1'-(Dc*C2+Cf)'  
       B1'*X+D21'*Bf'          B1'*R             -gam2*eye(nw)  -D21'*Dc'
       C1-Dc*C2                C1-(Dc*C2+Cf)     -Dc*D21        -eye(nz)       ];

end

