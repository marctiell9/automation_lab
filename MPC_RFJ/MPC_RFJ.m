%% MPC AUTOMATION LABORATORY
clear all
close all
clc

% we could consider in the cost function even other states, other than
% alpha absolute value
%% System model
% Model parameters
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
Vbar= 15; %maximum voltage

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
% C=[1 1 0 0];
C=[1 1 0 0;
    0 0 1 0]; % test weighting even theta_dot
D=zeros(1,1);

Ts  =   0.1;            % Sampling time   (increasing it improves performances)                    

Model           =   c2d(ss(A,B,C,D),Ts);  % Discrete-time model with zoh

% System matrices:
[Adt,Bdt,Cdt,Ddt]       =   ssdata(Model);                  % Model matrices 

%% Signal dimensions
nz     =   size(Adt,1);
nu      =   size(Bdt,2);
ny      =   size(Cdt,1);

%% Prediction horizon and cost function weighting matrices
N=5; % more demanding and worse response
Q       =   eye(ny);
Q(1,1)=30;
R       =   1e-4;

%% Inequality constraints
% Input inequalities
Au      =   [eye(nu);-eye(nu)];
bu      =   Vbar*ones(2*nu,1);
nqu     =   size(Au,1);                 % Number of input inequality constraints per stage

% State inequalities
Az     =   [1 0 0 0;-1 0 0 0];
bz     =   pi*ones(2,1); % -180<theta<180
nqz    =   size(Az,1);

%% Reference output and initial state
z0     =   zeros(nz,1);
% yref    =   pi/2*ones(ny,1);
yref = [pi/2*ones(1,1);
    zeros(1,1)];
%% Build overall matrices and vectors for QP (note - quadprog solves: min 0.5*x'*H*x + f'*x   subject to:  A*x <= b, Aeq*x = beq)
[Lambda_y,Gamma_y,Lambda_z,Gamma_z]   =   Traj_matrices(N,Adt,Bdt,Cdt,Ddt);
Qbar                                    =   zeros((N+1)*ny);
Rbar                                    =   zeros(N*nu);
Yref                                    =   zeros((N+1)*ny,1);                              
Aubar                                   =   zeros(N*nqu,N*nu);
bubar                                   =   zeros(N*nqu,1);
% Azbar                                  =   zeros((N+1)*nqz,(N+1)*nz); assenti perchè non inequality constraint sugli stati
% bzbar                                  =   zeros((N+1)*nqz,1);

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
Aeq     =   Gamma_z(end-nz+1:end,:)-Gamma_z(end-2*nz+1:end-nz,:);
beq     =   -(Lambda_z(end-nz+1:end,:)-Lambda_z(end-2*nz+1:end-nz,:))*z0;

% Cost function
f       =   z0'*Lambda_y'*Qbar*Gamma_y-Yref'*Qbar*Gamma_y;
H       =   (Gamma_y'*Qbar*Gamma_y)+Rbar;
H       =   0.5*(H+H');
                
%% QP options                   
options =   optimset('display','none');

%% Simulate with MPC
Tsimulation=10; %seconds [s]
Nsim               =   Tsimulation/Ts;
Zsim_MPC           =   zeros((Nsim+1)*nz,1);
Ysim_MPC           =   zeros(Nsim*ny,1);
Usim_MPC           =   zeros(Nsim*nu,1);
Zsim_MPC(1:nz,1)   =   z0;
zt                 =   z0;
tQP                =   zeros(Nsim-1,1);

U_prev=zeros(N,1);
% le matrici H e dei constraint rimangono costanti, variano le matrici
% delle traiettorie
for ind=2:Nsim+1
    bineq                               =   [bubar;bzbar-Azbar*Lambda_z*zt];
    beq                                 =   -(Lambda_z(end-nz+1:end,:)-Lambda_z(end-2*nz+1:end-nz,:))*zt;
    f                                   =   zt'*Lambda_y'*Qbar*Gamma_y-Yref'*Qbar*Gamma_y;
    tic
    U                                   =   quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],U_prev,options); % si puo usare il warm-starting con quad prog?
    tQP(ind-1,1)                        =   toc;
    U_prev=[U(2:N,1); U(N,1)];
    Usim_MPC((ind-2)*nu+1:(ind-1)*nu,1) =   U(1:nu,1);
    Zsim_MPC((ind-1)*nz+1:ind*nz,1)  =   Adt*Zsim_MPC((ind-2)*nz+1:(ind-1)*nz,1)+Bdt*Usim_MPC((ind-2)*nu+1:(ind-1)*nu,1);
    Ysim_MPC((ind-2)*ny+1:(ind-1)*ny,1) =   Cdt*Zsim_MPC((ind-2)*nz+1:(ind-1)*nz,1)+Ddt*Usim_MPC((ind-2)*nu+1:(ind-1)*nu,1);
    zt                                 =   Zsim_MPC((ind-1)*nz+1:ind*nz,1);
end
figure(),plot(Ts*[0:1:Nsim-1],Ysim_MPC(1:2:Nsim*ny)),title('alpha absolute')
figure(),plot(Ts*[0:1:Nsim],Zsim_MPC(1:4:Nsim*nz+4)),title('theta')
figure(),plot(Ts*[0:1:Nsim],Zsim_MPC(2:4:Nsim*nz+4)),title('alpha relative')
figure(),plot(Ts*[0:1:Nsim-1],Usim_MPC),title('Control input')
