clear
clc
%% MPC REAL TIME integral action
%% System model
Jl    =   3.28e-3;
Ks    =   1.3626;
Bl    =   0.0019;
Kg    =   70;
Beq   =   0.0072;   
Jeq   =   0.0023;
eta_m =   0.7029;
eta_g =   0.9012;
Kc    =   0.0083;
Kv    =   0.0085;
Rm    =   2.5241;
Vbar= 10; %maximum voltage

A=[0 0 1 0;
    0 0 0 1;
    0 Ks/Jeq -Beq/Jeq +Bl/Jeq;
    0 -Ks*(Jl+Jeq)/(Jl*Jeq) Beq/Jeq -Bl*(Jl+Jeq)/(Jl*Jeq)];
B=[0; 0 ; 1/Jeq; -1/Jeq];
A0=A;
B0=B;
B=eta_m*eta_g*Kc*Kg/Rm*B0;
A(3,3)=A0(3,3)-B0(3)*eta_m*eta_g*Kc*Kv*Kg^2/Rm;
A(4,3)=A0(4,3)-B0(4)*eta_m*eta_g*Kc*Kv*Kg^2/Rm;
C=[1 1 0 0];
D=zeros(1,1);

           % Sampling time  
Ts=0.1;

Model           =   c2d(ss(A,B,C,D),Ts);  % Discrete-time model with zoh

% System matrices:
[Adt,Bdt,Cdt,Ddt]       =   ssdata(Model);                  % Model matrices 

%% SYSTEM IN VELOCITY FORM
Aeq=[Adt zeros(4,1); -Cdt*Adt eye(1)];
Beq=[Bdt ; -C*Bdt];
Ceq=[0 0 0 0 1]; % the output is just the error, needed in the cost function (error on alpha absolute)
Deq=D;

%% Signal dimensions
nz     =   5;
nu      =   size(Beq,2);
ny      =   size(Ceq,1);

%% Prediction horizon and cost function weighting matrices
N       =   5;
Q       =   1;
R       =   1e-2;

%% Inequality constraints
% Input inequalities
Au      =   [eye(nu);-eye(nu)];
bu      =   10*ones(2*nu,1);
nqu     =   size(Au,1);                 % Number of input inequality constraints per stage

% State inequalities
Az     =   [1 0 0 0 0;-1 0 0 0 0];
bz     =   pi*ones(2,1); % -180<theta<180
nqz    =   size(Az,1);

%% Reference output and initial state
z0     =   zeros(nz,1);
yref    =   pi/2*ones(ny,1);

%% Build overall matrices and vectors for QP (note - quadprog solves: min 0.5*x'*H*x + f'*x   subject to:  A*x <= b, Aeq*x = beq)
[Lambda_y,Gamma_y,Lambda_z,Gamma_z]   =   Traj_matrices(N,Aeq,Beq,Ceq,Deq);
Qbar                                    =   zeros((N+1)*ny);
Rbar                                    =   zeros(N*nu);
Yref                                    =   zeros((N+1)*ny,1);                              
Aubar                                   =   zeros(N*nqu,N*nu);
bubar                                   =   zeros(N*nqu,1);
Azbar                                  =   zeros((N+1)*nqz,(N+1)*nz); 
bzbar                                  =   zeros((N+1)*nqz,1);

for ind = 1:N+1
    Qbar((ind-1)*ny+1:ind*ny,(ind-1)*ny+1:ind*ny)           =   Q;
    Yref((ind-1)*ny+1:ind*ny,1)                             =   yref;
    Azbar((ind-1)*nqz+1:ind*nqz,(ind-1)*nz+1:ind*nz)   =   Az;
    bzbar((ind-1)*nqz+1:ind*nqz,1)                       =   bz;
end

for ind = 1:N
    Rbar((ind-1)*nu+1:ind*nu,(ind-1)*nu+1:ind*nu)           =   R;
    Aubar((ind-1)*nqu+1:ind*nqu,(ind-1)*nu+1:ind*nu)        =   Au;
    bubar((ind-1)*nqu+1:ind*nqu,1)                          =   bu;
end

Aineq   =   [Aubar;Azbar*Gamma_z];
%Aineq   =   Aubar;
bineq   =   [bubar;bzbar-Azbar*Lambda_z*z0];
%bineq   =   bubar;

% Terminal equality constraint
% Aeq     =   Gamma_z(end-nz+1:end,:)-Gamma_z(end-2*nz+1:end-nz,:);
% beq     =   -(Lambda_z(end-nz+1:end,:)-Lambda_z(end-2*nz+1:end-nz,:))*z0;

% Cost function
f       =   z0'*Lambda_y'*Qbar*Gamma_y-Yref'*Qbar*Gamma_y;
H       =   (Gamma_y'*Qbar*Gamma_y)+Rbar;
H       =   0.5*(H+H');

%% QP options
options =   optimset('display','none');

%% Simulate with MPC
Tsimulation        =   10;         %seconds [s]
Nsim               =   Tsimulation/Ts;
Zsim_MPC           =   zeros((Nsim+1)*4,1);
Ysim_MPC           =   zeros(Nsim*ny,1);
deltaUsim_MPC           =   zeros(Nsim*nu,1);

z_temp=z0(1:4); %primo z_temp è proprio z0
z0_real=z0(1:4);
Zsim_MPC(1:4,1)   =   z0_real;
zt                 =   z0;
tQP                =   zeros(Nsim-1,1);
u=0;
uvec=zeros(Nsim-1,1);
deltaU_prev=zeros(N,1);

%% OBSERVER DESIGN
obs_poles=[-100,-110,-120,-130];
C=[1 0 0 0;0 1 0 0];
L=place(A',C',obs_poles);
eig(A-L'*C)
A_ob = A - L'*C;
B_ob = [B , L'];
C_ob = eye(4);




