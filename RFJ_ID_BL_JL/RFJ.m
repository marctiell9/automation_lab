function [zdot]=RFJ(t,z,u,th)
%% Read parameters, states and inputs
% Parameters
Jl       =       th(1,1);     %link inertia
Ks       =       th(2,1);     %for the moment it was supposed eqaul to zero, so it does not appear in the model
Bl      =       th(3,1);     % Torsional stiffness
Kg      =       th(4,1);
Beq      =       th(5,1);     
Jeq       =       th(6,1);     
eta_m      =     th(7,1);   % DC machine efficiency  
eta_g      =     th(8,1);    % gear box efficiency 
Kc      =       th(9,1);     % DC machine constants
Kv      =       th(10,1);    % 
Rm      =       th(11,1);    % DC machine resistance



%% Model equations


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
zdot= A*z + B*u;

