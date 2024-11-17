%% DATA RESET
clc
clear all

%possibile problema: i valori relativi al transitorio vanno scartati
%supponiamo di scartare i dati fino a 2 secondo
%% DATA LOADING

folder='data\19032024_SineII';
filePattern = fullfile(folder, '*.mat');
matFiles = dir(filePattern);
for k = 1:length(matFiles)
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(folder,baseFileName);
    fprintf(1,'%i ) Uploading %s \n',k,baseFileName);
    freqs(k)=str2double(extractBetween(baseFileName,'V','rads'));
    temp_file=load(fullFileName);
    data_names=fieldnames(temp_file);
    matData{k}=temp_file.(data_names{1});
end

fprintf('Data correctly uploaded \n');

%% SIGNAL FILTERING

fs=(1/0.002)*2*pi; %sampling frequency in radiants

for i = 1:length(matData)
    % conversion to radaints
    t{i}=matData{i}(1,:);
    v{i}=matData{i}(2,:);
    v{i}=-v{i};
    theta{i}=matData{i}(3,:);
    theta{i}=-theta{i}*2*pi/4096;
    alpha{i}=matData{i}(4,:);
    alpha{i}=alpha{i}*2*pi/4096;
    
    % we find the peak of theta and alpha and voltage
    [pks{i},locs{i}] = findpeaks(theta{i}); %locs vettore indice in cui ci sono picchi 
    [alpha_pks{i},alpha_locs{i}] = findpeaks(alpha{i});
    [v_pks{i},v_locs{i}] = findpeaks(v{i});
    

    n_pks=5;            % Number of peaks included in the transient, dopo questi picchi il transitorio viene cosnderato esaurito
    trans=locs{i}(n_pks)-round((locs{i}(n_pks)-locs{i}(n_pks-1))/4); %indidce del vettore nel quale finisce il transitotrio 
    mean_f=mean(theta{i}(1,locs{i}(n_pks-2):locs{i}(n_pks)));

    % ANTI-CAUSAL HIGHPASS BUTTERWORTH FILTER
    fpass=(freqs(i)-2.5)/(2*pi);
    [z,p,k] = butter(10,fpass*0.004,"high"); %high-pass filter
    sos = zp2sos(z,p,k);
    theta_f{i}=filtfilt(sos,1.02,theta{i});

    [pks_f{i},locs_f{i}]=findpeaks(theta_f{i},'MinPeakProminence',0.01); %picchi filtrati
    theta_mag(i)=mean(pks_f{i}(1,1+n_pks:end-n_pks));
    alpha_mag(i)=mean(alpha_pks{i}(1,1+n_pks:end-n_pks));
    v_mag(i)=mean(v_pks{i});

%     TO PLOT THETA WITH ITS PEAKS
%     
    % figure(i)
    % plot(t{i},alpha{i},t{i}(alpha_locs{i}),alpha_pks{i},"o")
    % xlabel("time(sec)")
    % ylabel("alpha(rad)")
    % axis tight

    mag_h_theta(i)=theta_mag(i)/v_mag(i);
    mag_h_alpha(i)=alpha_mag(i)/v_mag(i);

    phase_h_theta(i)=-180+phdiffmeasure(v{i}(1,locs{i}(n_pks):locs{i}(end-n_pks)),theta_f{i}(1,locs{i}(n_pks):locs{i}(end-n_pks)),fs/(2*pi),'dft')*180/pi;
    phase_h_alpha(i)=180+phdiffmeasure(v{i}(1,locs{i}(n_pks):locs{i}(end-n_pks)),alpha{i}(1,locs{i}(n_pks):locs{i}(end-n_pks)),fs/(2*pi),'dft')*180/pi;

    theta_f{i}=[theta{i}(1,1:(trans-1)),theta_f{i}(1,trans:end)+mean_f]; %il segnale filtrato finale è questo, perchè durante il transitorio e dopo viene trattato diversamente
    %TO PLOT ORIGINAL THETA VS FILTERED THETA
    % figure(i)
    % plot(t{i},theta{i},'LineWidth',1.5);
    % hold on
    % plot(t{i},theta_f{i},'r','LineWidth',1.5);
    % legend('Original','Detrended','Location','southeast');
    % hold on

end

theta_mag(6)=0.35;      % Correction on lowest freq. magnitude, messo a mano perchè non funziona in maniera automatica

% mag_h_theta=mag2db(mag_h_theta);
% mag_h_alpha=mag2db(mag_h_alpha);

% hold off;



%% SYSTEM MODEL 

%  DC Motor Parameters
R_m=2.5241;      % Motor armature resistence
eta_g=0.9012;    % Gearbox Efficiency
eta_m=0.7029;    % Motor Efficiency
K_g=70;          % High-gear Total Gear Ratio
k_t=0.0083;  % Motor Current-Torque Constant
k_m=0.0085;  % Motor Back-EMF Constant

%  Gears Parameters
m24=0.005;       % 24-thooth gear mass
m120=0.083;      % 120-thooth gear mass
r24=6.35*10^-3;  % 24-thooth gear radius
r120=0.032;      % 120-thooth gear radius

%  Rotoflex Parameters
L1=29.8*1e-2;                                % Main Arm Lenght
L2=15.6*1e-2;                                % Load Arm Lenght
m1=0.064;                                    % Main Arm Mass
m2=0.03;                                     % Load Arm Mass
d12=21*1e-2;                                 % Joint to middle of Load Arm
J_l=(m1*(L1^2)/3)+(m2*(L2^2)/12)+(m2*d12^2); % Equivalent Moment of Inertia of Arm
B_l=0.00316;                                 % Equivalent Viscous Damping Coefficient of the Arm
B_l=0.0019;
J_eq=0.0023;                                 % High-Gear Equivalent Moment of Inertia with No Load
B_eq=0.0072;                                 % High-Gear Equivalent Viscous Damping Coefficient with No Load
K_s=1.3597;                                  % Torsional Spring Stiffness (Taken from complete report)
K_s=1.3626;

%  Rotoflex Linear State-Space Model
A=[0                0               1                   0;
   0                0               0                   1;
   0            (K_s/J_eq)    -(B_eq/J_eq)          (B_l/J_eq);
   0 -K_s*(J_l+J_eq)/(J_l*J_eq) B_eq/J_eq -B_l*(J_l+J_eq)/(J_l*J_eq)];
B=[0;0;1/J_eq;-1/J_eq];
C=[1 0 0 0;
    0 1 0 0];
D=0;
A0=A;
B0=B;
B=eta_m*eta_g*k_t*K_g/R_m*B0;
A(3,3)=A0(3,3)-B0(3)*eta_m*eta_g*k_t*k_m*K_g^2/R_m;
A(4,3)=A0(4,3)-B0(4)*eta_m*eta_g*k_t*k_m*K_g^2/R_m;
sys=ss(A,B,C,D);
%% BODE PLOTS

% bode(sys,freqs,'.-')
% grid on

w = logspace(0,3,2000);

[mag,phase] = bode(sys,w);
mag_theta=mag(1,1,:);
mag_alpha=mag(2,1,:);
phase_theta=phase(1,1,:);
phase_alpha=phase(2,1,:);

figure(1)
subplot(4,1,1)
loglog(w, squeeze(mag_theta), freqs, mag_h_theta, 'ro','LineWidth',1.5)
xlabel('freq(rad/s)');
ylabel('mag[\theta/V]');
grid
subplot(4,1,2)
semilogx(w, squeeze(phase_theta), freqs, phase_h_theta, 'ro','LineWidth',1.5)
xlabel('freq(rad/s)');
ylabel('phase[\theta/V](deg)');
grid
subplot(4,1,3)
loglog(w, squeeze(mag_alpha), freqs, mag_h_alpha, 'ro','LineWidth',1.5)
xlabel('freq(rad/s)');
ylabel('mag[\alpha/V]');
grid
subplot(4,1,4)
semilogx(w, squeeze(phase_alpha), freqs,[phase_h_alpha(1,1:6) phase_h_alpha(1,7)-360 phase_h_alpha(1,8:end)], 'ro','LineWidth',1.5)
xlabel('freq(rad/s)');
ylabel('phase[\alpha/V](deg)');
grid


figure()
loglog(w, squeeze(mag_theta), freqs, mag_h_theta, 'ro','LineWidth',1.5)
xlabel('freq(rad/s)');
ylabel('mag[\theta/V]');

figure()
semilogx(w, squeeze(phase_alpha), freqs,[phase_h_alpha(1,1:6) phase_h_alpha(1,7)-360 phase_h_alpha(1,8:end)], 'ro','LineWidth',1.5)
xlabel('freq(rad/s)');
ylabel('mag[\theta/V]');
%% MIDA2 
% theta_f1=cell2mat(theta_f);
% N=length(theta_f{1}); 
% Ts=0.002; 
% Tend=(N-1)*Ts;
% time=0:Ts:Tend;
% 
% 
% 
% w=freqs; 
% 
% for i=1:1:length(w)
% %Matrix A is unique both to retrive alpha and theta frequency response
% A=[ones(1,N)*(sin(w(i).*time).*sin(w(i).*time))'  ones(1,N)*(sin(w(i).*time).*cos(w(i).*time))';
%    ones(1,N)*(sin(w(i).*time).*cos(w(i).*time))'  ones(1,N)*(cos(w(i).*time).*cos(w(i).*time))' ];
% 
% % Instead matrix b is different beacuse it contains the measurment of alpha and theta
% 
% b_theta=[ones(1,N)*(theta_f{i}.*sin(w(i).*time))';
%         ones(1,N)*(theta_f{i}.*cos(w(i).*time))'];
% 
% b_alpha=[ones(1,N)*(alpha{i}.*sin(w(i).*time))';
%         ones(1,N)*(alpha{i}.*cos(w(i).*time))'];
% 
% temp_theta=A^-1*b_theta;
% a_hat_theta(i)=temp_theta(1);
% b_hat_theta(i)=temp_theta(2);
% 
% temp_alpha=A^-1*b_alpha;
% a_hat_alpha(i)=temp_alpha(1);
% b_hat_alpha(i)=temp_alpha(2);
% 
% phi_theta(i)=atan(b_hat_theta(i)/a_hat_theta(i));
% phi_alpha(i)=atan(b_hat_alpha(i)/a_hat_alpha(i));
% B_theta(i)=(a_hat_theta(i)/cos(phi_theta(i))+b_hat_theta(i)/sin(phi_theta(i)))/2;
% B_alpha(i)=(a_hat_alpha(i)/cos(phi_alpha(i))+b_hat_alpha(i)/sin(phi_alpha(i)))/2;
% 
% end
% C_theta=B_theta/v_mag;
% C_alpha=B_alpha/v_mag;
% 
% 
% w = logspace(0,3,2000);
% 
% [mag,phase] = bode(sys,w);
% mag_theta=mag(1,1,:);
% mag_alpha=mag(2,1,:);
% phase_theta=phase(1,1,:);
% phase_alpha=phase(2,1,:);
% 
% figure(1)
% subplot(4,1,1)
% loglog(w, squeeze(mag_theta), freqs, C_theta, 'ro','LineWidth',1.5)
% xlabel('freq(rad/s)');
% ylabel('mag[\theta/V]');
% grid
% subplot(4,1,2)
% semilogx(w, squeeze(phase_theta), freqs, rad2deg(phi_theta)-180, 'ro','LineWidth',1.5)
% xlabel('freq(rad/s)');
% ylabel('phase[\theta/V](deg)');
% grid
% subplot(4,1,3)
% loglog(w, squeeze(mag_alpha), freqs, C_alpha, 'ro','LineWidth',1.5)
% xlabel('freq(rad/s)');
% ylabel('mag[\alpha/V]');
% grid
% subplot(4,1,4)
% semilogx(w, squeeze(phase_alpha), freqs, rad2deg(phi_alpha)+180, 'ro','LineWidth',1.5)
% xlabel('freq(rad/s)');
% ylabel('phase[\alpha/V](deg)');
% grid