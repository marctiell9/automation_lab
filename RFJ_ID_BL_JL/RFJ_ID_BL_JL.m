clear all
close all
clc
%% Model parameters
%Now we suppose to know Jl, Ks, and Bl
Jl       =  3.28e-3;     % Link inertia
Bl=0.0019;
Ks=1.3626;
Kg=70;
% Beq      =       th(4,1);     
% Jeq       =       th(5,1);     
% eta_m      =     th(6,1);   % DC machine efficiency  
% eta_g      =     th(7,1);    % gear box efficiency 
% Kc      =       th(8,1);     % DC machine constants
% Kv      =       th(9,1);    % 
% Rm      =       th(10,1);    % DC machine resistance
th=[Jl;Ks;Bl;Kg];
%% Load data
% %load 2.5Volt_1Hertz_try02.mat
% load 5.0Volt_1Hertz_try01.mat
% % load 0-2Hz_sweep_1V.mat
% load 1VPulse_1Sec.mat
% %The sampling time used for measuring was 0.002
% downsampling 	=   1;                                                     % Reduce sampling frequency for higher simulation speed with FFD
% Ts             	=   downsampling*0.002;   
% Tin=1;
% 
% % Square wave 5V
% uin             =   data_01_Mar_2024_12_31_23(2,1:downsampling:end).*-1; %-1 factor HERE, it should be wrong   
% ymeas          	=   [data_01_Mar_2024_12_31_23(3,1:downsampling:end).*2*pi/4096';   % the factor 360/4096 is needed to go form the encoder measure into radiants(or degree?)     
%                  	 data_01_Mar_2024_12_31_23(4,1:downsampling:end).*-1*2*pi/4096'];   %-1 factor here to flap alpha measures
% 
% %1V PULSE
% Tin=500;
% Tend=8000;
% uin             =   data_12_Mar_2024_15_36_23(2,1:downsampling:Tend).*-1;    
% ymeas          	=   [data_12_Mar_2024_15_36_23(3,1:downsampling:Tend).*2*pi/4096';   % the factor 360/4096 is needed to go form the encoder measure into radiants(or degree?)     
%                    	 data_12_Mar_2024_15_36_23(4,1:downsampling:Tend).*-2*pi/4096'];
% 
% z0=[ymeas(1,Tin);0;ymeas(2,Tin);0];
% 
% %1V PULSE
% Tin=500;
% Tend=8000;
% uin             =   data_12_Mar_2024_15_36_23(2,1:downsampling:Tend).*-1;    
% ymeas          	=   [data_12_Mar_2024_15_36_23(3,1:downsampling:Tend).*2*pi/4096';   % the factor 360/4096 is needed to go form the encoder measure into radiants(or degree?)     
%                    	 data_12_Mar_2024_15_36_23(4,1:downsampling:Tend).*-2*pi/4096'];
% 
% z0=[ymeas(1,Tin);0;ymeas(2,Tin);0];

%% NEW DATA
%Square wave 2V 5rad/s
downsampling 	=   1;                                                     % Reduce sampling frequency for higher simulation speed with FFD
Ts             	=   downsampling*0.002;   
Tin=1;
%data=load('2V5rads.mat') ;
data=load('2.5Volt_1Hertz_try02.mat');
data_names =  fieldnames(data);
uin             =   data.(data_names{1})(2,1:downsampling:end).*-1;    
ymeas          	=   [data.(data_names{1})(3,1:downsampling:end).*2*pi/4096';      
                   	 data.(data_names{1})(4,1:downsampling:end).*-1*2*pi/4096'];
figure,plot(data.(data_names{1})(1,1:downsampling:end),uin(1,:),'LineWidth',2),xlabel('time [s]'),ylabel('Volt [V]'),title('Input signal'),grid
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:).*180/pi),xlabel('time [s]'),ylabel('degree[°]'),title('Theta experimental data '),grid
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:).*180/pi),xlabel('time [s]'),ylabel('degree[°]'),title('Alpha experimental data '),grid


% [pks{i},locs{i}] = findpeaks(theta{i}); %locs vettore indice in cui ci sono picchi 
% [alpha_pks{i},alpha_locs{i}] = findpeaks(alpha{i});
% [v_pks{i},v_locs{i}] = findpeaks(v{i});
% n_pks=5;
% trans=locs{i}(n_pks)-round((locs{i}(n_pks)-locs{i}(n_pks-1))/4);
% mean_f=mean(theta{i}(1,locs{i}(n_pks-2):locs{i}(n_pks)));
% % ANTI-CAUSAL HIGHPASS BUTTERWORTH FILTER
% fpass=(freqs(i)-2.5)/(2*pi);
% [z,p,k] = butter(10,fpass*0.004,"high"); %high-pass filter
% sos = zp2sos(z,p,k);
% theta_f{i}=filtfilt(sos,1.02,theta{i});
% [pks_f{i},locs_f{i}]=findpeaks(theta_f{i},'MinPeakProminence',0.01); %picchi filtrati
% theta_f{i}=[theta{i}(1,1:(trans-1)),theta_f{i}(1,trans:end)+mean_f];

z0=[ymeas(1,Tin);0;ymeas(2,Tin);0];

%% Select weight matrix for cost function and scaling factors
Q               =  1*eye(2);
Q(2,2)=10;
%Q(1,1)=100;
scaling         =   ones(7,1);
%scaling =[1e2;1e3;1e1;1e1;1e3;1e3;1e-1;1];
%% Initialize parameter estimate
x0             	=   [0.015;2.08e-3;0.69;0.9;7.68e-3;7.68e-3;2.6].*scaling;                      % intial values of the optimization variables will be the ones of the data sheet
%x0 = zeros(10,1);
%% Constraints
% we have 10 parameters to estimate
% Maybe Kg (gear ratio can be assumed to be known)
C       =       [ 1 zeros(1,6);                              %Beq>0                           
                  0  1 zeros(1,5);                          % Jeq
                  0  -1 zeros(1,5);                         %Jeq
                   zeros(1,2) 1 zeros(1,4);                   %eta_m
                   zeros(1,2) -1 zeros(1,4);
                   zeros(1,3) 1 zeros(1,3);                    %eta_g
                   zeros(1,3) -1 zeros(1,3);
                   zeros(1,4) 1 zeros(1,2);                    %Kc
                   zeros(1,4) -1 zeros(1,2);
                   zeros(1,5) 1 zeros(1,1);                    %Kv
                   zeros(1,5) -1 zeros(1,1);
                   zeros(1,6) 1;                              %Rm
                   zeros(1,6) -1; 
                    ];
d       =       [0; 2.087e-3-2.087e-3*10/100; -2.087e-3-2.087e-3*10/100;0.69-0.69*5/100; -0.69-0.69*5/100; 0.9-0.9/100 ; -0.9-0.9/100 ; 7.68e-3-7.68e-3*12/100 ; -7.68e-3-7.68e-3*12/100 ; 7.68e-3-7.68e-3*12/100 ; -7.68e-3-7.68e-3*12/100 ; 2.6-2.6*12/100;-2.6-12*10/100 ];


%% Gauss-Newton
% Initialize solver options
myoptions               =   myoptimset;
myoptions.Hessmethod  	=	'GN'; 
myoptions.GN_funF       =	@(x)RFJ_sim_err(x,z0,uin,ymeas,th,Ts,Q,scaling); 
myoptions.gradmethod  	=	'CD';
myoptions.graddx        =	2^-17;
myoptions.tolgrad    	=	1e-12;
myoptions.GN_sigma      =	0;
myoptions.ls_nitermax   =	1e3;
myoptions.nitermax      =	1000;
myoptions.tolfun        =   1e-18;
myoptions.ls_beta    	=	0.9;    
%myoptions.ls_c          =	1;
outputfun   =   myoptions.outputfcn;
%% Run solver
%[xstar,fxstar,niter,exitflag,xsequence] = myfmincon(@(x)RFJ_sim_cost(x,z0,uin,ymeas,th,Ts,Q,scaling),x0,[],[],C,d,0,0,myoptions);
%% Code generation
codegen RFJ_sim_cost -args {x0,z0,uin,ymeas,th,Ts,Q,scaling}
 [xstar,fxstar,niter,exitflag,xsequence] = myfmincon(@(x)RFJ_sim_cost_mex(x,z0,uin,ymeas,th,Ts,Q,scaling),x0,[],[],C,d,0,0,myoptions);
%% PULSE IDENTIFICATION
% [xstar./scaling]
% [~,ysim]               	=   RFJ_sim_cost(xstar,z0,uin,ymeas,th,Ts,Q,scaling);
% figure,plot(data_12_Mar_2024_15_36_23(1,1:downsampling:Tend),ymeas(1,:),data_12_Mar_2024_15_36_23(1,1:downsampling:Tend),ysim(1,:)),legend('measure','ID model'),title('theta')
% figure,plot(data_12_Mar_2024_15_36_23(1,1:downsampling:Tend),ymeas(2,:),data_12_Mar_2024_15_36_23(1,1:downsampling:Tend),ysim(2,:)),legend('measure','ID model'),title('alpha')

%% SQUARE WAVE 5V 1hertz identification
[xstar./scaling]
[~,ysim]               	=   RFJ_sim_cost(xstar,z0,uin,ymeas,th,Ts,Q,scaling);
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:).*180/pi,data.(data_names{1})(1,1:downsampling:end),ysim(1,:).*180/pi),legend('Measured theta','Modeled theta' ),title('theta'),xlabel('time [s]'),ylabel('degree[°]'),title('Modeled theta vs Measured theta'),grid
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:).*180/pi,data.(data_names{1})(1,1:downsampling:end),ysim(2,:).*180/pi),legend('Measured alpha','Modeled alpha'),title('alpha'),xlabel('time [s]'),ylabel('degree[°]'),title('Modeled alpha vs Measured alpha'),grid

%% VALIDATION PULSE'Measured alpha',
% data=load('1VPulse_1Sec.mat'); %peggiore
 data=load('2.5VPulse_0.5sec.mat'); %meglio
%data=load('5V1s.mat'); %forse meglio meglio
downsampling=1;
data_names =  fieldnames(data);
uin             =   data.(data_names{1})(2,1:downsampling:end).*-1;    
ymeas          	=   [data.(data_names{1})(3,1:downsampling:end).*2*pi/4096';      
                   	 data.(data_names{1})(4,1:downsampling:end).*-1*2*pi/4096'];
figure,plot(data.(data_names{1})(1,1:downsampling:end-1),uin(1,1:end-1),'LineWidth',3),xlabel('time [s]'),ylabel('volt [V]'),title('Input signal')
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:).*180/pi),xlabel('time [s]'),ylabel('degree [°]'),title('Validation data theta')
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:).*180/pi),xlabel('time [s]'),ylabel('degree [°]'),title('Validation data alpha')
z0=[ymeas(1,1);0;ymeas(2,1);0];
%th=[Jl;Ks;Bl;Kg;xstar./scaling];
%ultimi parametri trovati
th=[Jl;Ks;Bl;Kg;Beq;Jeq;eta_g;eta_m;Kc;Kv;Rm];
Jl       =  3.28e-3;     % Link inertia
Bl=0.0019;
Ks=1.3626;
Kg=70;
Beq=0.0041;
Jeq=   0.0023;
eta_g=  0.6916;
eta_m=0.9003;
Kc=0.0085;
Kv=0.0079;
Rm=2.9077;

N_FFD               =       size(ymeas,2);      % number of samples
zhat            =       zeros(4,N_FFD);     
zhat(:,1)       =       z0;

for ind=2:N_FFD
    zdot               =   RFJ(0,zhat(:,ind-1),uin(1,ind-1),th);
    zhat(:,ind)    =   zhat(:,ind-1)+Ts*zdot;
end

figure,plot(data.(data_names{1})(1,1:downsampling:end),uin(1,:)),title('pulse validation')
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:).*180/pi,data.(data_names{1})(1,1:downsampling:end),zhat(1,:).*180/pi),legend('Measured theta','Modeled theta'),xlabel('time [s]'),ylabel('degree[°]'),title('Modeled theta vs Measured theta')
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:).*180/pi,data.(data_names{1})(1,1:downsampling:end),zhat(2,:).*180/pi),legend('Measured alpha','Modeled alpha'),xlabel('time [s]'),ylabel('degree[°]'),title('Modeled alpha vs Measured alpha')
% figure,plot(data_01_Mar_2024_12_26_39(1,1:downsampling:end),ymeas(1,:)-zhat(1,:)),title('theta error')
% figure,plot(data_01_Mar_2024_12_26_39(1,1:downsampling:end),ymeas(2,:)-zhat(2,:)),title('alpha error')

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
C=[1 0 0 0;
    0 1 0 0];
D=zeros(2,1);
x0=zeros(1,4);
sys_id=ss(A,B,C,D); 
s=tf('s');
G=C*(s*eye(4)-A)^-1*B; 
Gs=minreal(G,1e-6);
Gtheta=Gs(1);
theta_sim=lsim(sys_id(1),uin,data.(data_names{1})(1,1:downsampling:end));
alpha_sim=lsim(sys_id(2),uin,data.(data_names{1})(1,1:downsampling:end));
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:),data.(data_names{1})(1,1:downsampling:end),theta_sim),legend('Theta measure','theta simulation')
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:),data.(data_names{1})(1,1:downsampling:end),alpha_sim),legend('Theta measure','alpha simulation')



%% VALIDATION Square wave 2V 5rad/s (problem with shift)
% data=load('1VPulse_1Sec.mat'); %peggiore
% data=load('2.5VPulse_0.5sec.mat'); %meglio
data=load('2V5rads.mat'); %forse meglio meglio
data_names =  fieldnames(data);
uin             =   data.(data_names{1})(2,1:downsampling:end).*-1;    
ymeas          	=   [data.(data_names{1})(3,1:downsampling:end).*2*pi/4096';      
                   	 data.(data_names{1})(4,1:downsampling:end).*-1*2*pi/4096'];
figure,plot(data.(data_names{1})(1,1:downsampling:end),uin(1,:))
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:))
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:))
z0=[ymeas(1,1);0;ymeas(2,1);0];
th=[Jl;Ks;Bl;Kg;xstar./scaling];
N_FFD               =       size(ymeas,2);      % number of samples
zhat            =       zeros(4,N_FFD);     
zhat(:,1)       =       z0;

for ind=2:N_FFD
    zdot               =   RFJ(0,zhat(:,ind-1),uin(1,ind-1),th);
    zhat(:,ind)    =   zhat(:,ind-1)+Ts*zdot;
end

figure,plot(data.(data_names{1})(1,1:downsampling:end),uin(1,:)),title('Square wave 2V 5rad/s validation')
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:),data.(data_names{1})(1,1:downsampling:end),zhat(1,:)),legend('measure','ID model'),title('theta')
figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:),data.(data_names{1})(1,1:downsampling:end),zhat(2,:)),legend('measure','ID model'),title('alpha')




%% RESULTS SINE WAVE 
% [xstar./scaling]
% [~,ysim]               	=   RFJ_sim_cost(xstar,z0,uin,ymeas,th,Ts,Q,scaling);
% figure,plot(data_12_Mar_2024_15_56_32(1,Tin:downsampling:end),ymeas(1,:),data_12_Mar_2024_15_56_32(1,Tin:downsampling:end),ysim(1,:)),legend('measure','ID model'),title('theta')
% figure,plot(data_12_Mar_2024_15_56_32(1,Tin:downsampling:end),ymeas(2,:),data_12_Mar_2024_15_56_32(1,Tin:downsampling:end),ysim(2,:)),legend('measure','ID model'),title('alpha')



%% VALIDATION

% SINE WAVE
% load 2V_3Hz_sine.mat
% uin             =   data_12_Mar_2024_15_56_32(2,1:downsampling:end).*-1;    
% ymeas          	=   [data_12_Mar_2024_15_56_32(3,1:downsampling:end).*2*pi/4096';      
%                    	 data_12_Mar_2024_15_56_32(4,1:downsampling:end).*-1*2*pi/4096'];
% theta_mean=mean(ymeas(1,1:downsampling:end));
% alpha_mean=mean(ymeas(2,1:downsampling:end));
% ymeas=ymeas-[theta_mean; alpha_mean];
% order = 1;              % Filter order
% cutoffFreq = 1;       % Cutoff frequency in Hz
% samplingFreq = 500;    % Sampling frequency in Hz
% [b, a] = butter(order, cutoffFreq/(samplingFreq/2), 'high');
% ymeas(1,:)=filtfilt(b, a, ymeas(1,:));
% z0=[ymeas(1,1);0;ymeas(2,1);0];
% th=[Jl;Ks;Bl;Kg;xstar./scaling];

% N_FFD               =       size(ymeas,2);      % number of samples
% zhat            =       zeros(4,N_FFD);     
% zhat(:,1)       =       z0;
% 
% for ind=2:N_FFD
%     zdot               =   RFJ(0,zhat(:,ind-1),uin(1,ind-1),th);
%     zhat(:,ind)    =   zhat(:,ind-1)+Ts*zdot;
% end
% 
% theta_mean_hat=mean(zhat(1,1:downsampling:end));
% alpha_mean_hat=mean(zhat(2,1:downsampling:end));
% zhat(1:2,:) = zhat(1:2,:)-[theta_mean_hat; alpha_mean_hat];
% 
% figure,plot(data_12_Mar_2024_15_56_32(1,1:downsampling:end),uin(1,:)),title('sine wave validation')
% figure,plot(data_12_Mar_2024_15_56_32(1,1:downsampling:end),ymeas(1,:),data_12_Mar_2024_15_56_32(1,1:downsampling:end),zhat(1,:)),legend('measure','ID model'),title('theta')
% figure,plot(data_12_Mar_2024_15_56_32(1,1:downsampling:end),ymeas(2,:),data_12_Mar_2024_15_56_32(1,1:downsampling:end),zhat(2,:)),legend('measure','ID model'),title('alpha')
% figure,plot(data_12_Mar_2024_15_56_32(1,1:downsampling:end),ymeas(1,:)-zhat(1,:)),title('theta error')
% figure,plot(data_12_Mar_2024_15_56_32(1,1:downsampling:end),ymeas(2,:)-zhat(2,:)),title('alpha error')
% 

%% Validation square wave
% uin             =   data_01_Mar_2024_12_26_39(2,1:downsampling:end).*-1;    
% ymeas          	=   [data_01_Mar_2024_12_26_39(3,1:downsampling:end).*2*pi/4096';      
%                    	 data_01_Mar_2024_12_26_39(4,1:downsampling:end).*-1*2*pi/4096'];
% theta_mean=mean(ymeas(1,1:downsampling:end));
% alpha_mean=mean(ymeas(2,1:downsampling:end));
% ymeas=ymeas-[theta_mean; alpha_mean];
% z0=[ymeas(1,1);0;ymeas(2,1);0];
% th=[Jl;Ks;Bl;Kg;xstar./scaling];
% N_FFD               =       size(ymeas,2);      % number of samples
% zhat            =       zeros(4,N_FFD);     
% zhat(:,1)       =       z0;
% 
% for ind=2:N_FFD
%     zdot               =   RFJ(0,zhat(:,ind-1),uin(1,ind-1),th);
%     zhat(:,ind)    =   zhat(:,ind-1)+Ts*zdot;
% end
% 
% figure,plot(data_01_Mar_2024_12_26_39(1,1:downsampling:end),uin(1,:)),title('square wave validation')
% figure,plot(data_01_Mar_2024_12_26_39(1,1:downsampling:end),ymeas(1,:),data_01_Mar_2024_12_26_39(1,1:downsampling:end),zhat(1,:)),legend('measure','ID model'),title('theta')
% figure,plot(data_01_Mar_2024_12_26_39(1,1:downsampling:end),ymeas(2,:),data_01_Mar_2024_12_26_39(1,1:downsampling:end),zhat(2,:)),legend('measure','ID model'),title('alpha')
% figure,plot(data_01_Mar_2024_12_26_39(1,1:downsampling:end),ymeas(1,:)-zhat(1,:)),title('theta error')
% figure,plot(data_01_Mar_2024_12_26_39(1,1:downsampling:end),ymeas(2,:)-zhat(2,:)),title('alpha error')
% 
% 
% 


%% Validation sine sweep
% data= load('0-2Hz_sweep_1V.mat');
% 
% data_names =  fieldnames(data);
% 
% 
% 
% uin             =   data.(data_names{1})(2,1:downsampling:end).*-1;    
% ymeas          	=   [data.(data_names{1})(3,1:downsampling:end).*2*pi/4096';      
%                    	 data.(data_names{1})(4,1:downsampling:end).*-1*2*pi/4096'];
% theta_mean=mean(ymeas(1,1:downsampling:end));
% alpha_mean=mean(ymeas(2,1:downsampling:end));
% ymeas=ymeas-[theta_mean; alpha_mean];
% z0=[ymeas(1,1);0;ymeas(2,1);0];
% th=[Jl;Ks;Bl;Kg;xstar./scaling];
% N_FFD               =       size(ymeas,2);      % number of samples
% zhat            =       zeros(4,N_FFD);     
% zhat(:,1)       =       z0;
% 
% for ind=2:N_FFD
%     zdot               =   RFJ(0,zhat(:,ind-1),uin(1,ind-1),th);
%     zhat(:,ind)    =   zhat(:,ind-1)+Ts*zdot;
% end
% 
% theta_mean_hat=mean(zhat(1,1:downsampling:end));
% alpha_mean_hat=mean(zhat(2,1:downsampling:end));
% zhat(1:2,:) = zhat(1:2,:)-[theta_mean_hat; alpha_mean_hat];
% 
% figure,plot(data.(data_names{1})(1,1:downsampling:end),uin(1,:)),title('sine sweep validation')
% figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:),data.(data_names{1})(1,1:downsampling:end),zhat(1,:)),legend('measure','ID model'),title('theta')
% figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:),data.(data_names{1})(1,1:downsampling:end),zhat(2,:)),legend('measure','ID model'),title('alpha')
% figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(1,:)-zhat(1,:)),title('theta error')
% figure,plot(data.(data_names{1})(1,1:downsampling:end),ymeas(2,:)-zhat(2,:)),title('alpha error')




