function    [err_vec,ysim]    =   RFJ_sim_err(x,z0,uin,ymeas,th,Ts,Q,scaling)
% Function that computes the output trajectory of the RFJ model
% introduced in Lab. session A and the weighted error between the simulated output
% and the measured one. Variable x corresponds to the vehicle parameters 
% to be estimated, possibly scaled. The simulation is carried out with forward finite
% differences

%% Model parameters
x       =       x./scaling;
Jl       =       th(1,1);  
Ks       =       th(2,1);     %for the moment it was supposed eqaul to zero, so it does not appear in the model
Bl     =        th(3,1);     % Torsional stiffness
Kg=            th(4,1);    % Gearbox ratio
Beq       =      x(1,1);     % Motor inertia
Jeq      =       x(2,1);     
eta_m      =     x(3,1);   % DC machine efficiency  
eta_g      =     x(4,1);    % gear box efficiency 
Kc      =       x(5,1);     % DC machine constants
Kv      =       x(6,1);    % 
Rm      =       x(7,1);    % DC machine resistance
th      =       [Jl;Ks;Bl;Kg;Beq;Jeq;eta_m;eta_g;Kc;Kv;Rm];

%% Initialize simulation output
N               =       size(ymeas,2);      % number of samples
zhat            =       zeros(4,N);     % matrix with states
zhat(:,1)       =       z0;

% Run simulation with FFD
% for ind=2:N
%     zdot               =   RFJ(0,zhat(:,ind-1),uin(:,ind-1),th);
%     zhat(:,ind)    =   zhat(:,ind-1)+Ts*zdot;
% end

%% ODE45
% ODE45
for ind=2:N
    [t,zhat_temp]          =   ode45(@(t,z)RFJ(t,z,uin(:,ind-1),th),[0 Ts],zhat(:,ind-1));
    zhat(:,ind)            =   zhat_temp(end,:)';
end


%% Collect simulated output 
ysim    =   [zhat(1,:);zhat(2,:)]; %The outputs are theta and alpha only

%% Compute weighted errors
err     =   Q*(ymeas-ysim);

%% Stack errors in one vector
err_vec   =   [err(1,:) err(2,:)]'/sqrt(N); % squared error case (L2 NORM)
% err_vec =   [err(1,:) err(2,:)]'/sqrt(N)'; % l-inf norm case
