clc

SYS.A{1}=[0 -1-0.3 ; 1 -0.5];
SYS.B1{1}=[-2 0 ; 1 0];
SYS.C1{1}=[1 0];
SYS.C2{1}=100*[-1 1];
SYS.D21{1}=[0 1];
SYS.A{2}=[0 -1+0.3 ; 1 -0.5];
SYS.B1{2}=SYS.B1{1};
SYS.C1{2}=SYS.C1{1};
SYS.C2{2}=SYS.C2{1};
SYS.D21{2}=SYS.D21{1};

[FILTER,gamma,feas] = Hinf_RobustFilter( SYS )
