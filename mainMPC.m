%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This software performs battery charge control using Split-future MPC
%
% Copyright (c) 2021 by Aloisio Kawakita de Souza of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a MIT license. It is provided "as is", without express or implied
% warranty, for educational and informational purposes only.
%
% This file is provided as a supplement to: 
%
% M. A. Xavier, A. K. de Souza, and M. S. Trimboli, “A split-future MPC algorithm
% for lithium-ion battery cell-level fast-charge control,” in Preprints of the 21st IFAC
% World Congress (Virtual), pp. 12638– 12643, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


clear;close all;
tic
addpath('./functions');


% Load battery model for a Ford C-max battery (Panasonic 25 Ah NMC Prismatic cell)
load('FORD53model.mat');


Nsim = 800;      % Duration of simulation
time = 1:Nsim+1;  % Time vector

SOC_ref    = 90;      % SOC reference [%]
deltaT = 1;           % Sampling period

% Plant Initical Conditions
SOC0_plant = 0.1;           % Plant's real SOC
ir0_plant = 0;
h0_plant = 0;

% MPC Tuning variables
mpcData.flag = 'charge';
mpcData.Np = 20;    % Prediction horizon
mpcData.Nc = 2;     % Control horizon
mpcData.Ru  = 1e-6*eye(mpcData.Nc,mpcData.Nc);  % Input weighting
mpcData.Sigma = tril(ones(mpcData.Nc,mpcData.Nc));

u0 = 0;                     % Initial current
mpcData.uk_1 = u0;           % u[k-1]
mpcData.SOCk_1 = 0;         % SOC[k-1]
mpcData.DUk_1 = 0;
mpcData.lambda = [];

% MPC constraints
mpcData.const.u_max = 200;    % Max discharge current
mpcData.const.u_min = -150;   % Max charge current
mpcData.const.v_min = 3.2;    % Min voltage
mpcData.const.v_max = 4.2;    % Max voltage
mpcData.const.z_max = SOC_ref/100;  % SOC reference during charging
mpcData.const.z_min = SOC_ref/100;   % SOC reference during discharging
mpcData.model = model;       % Store model to use later
mpcData.deltaT = deltaT;

% Initialize plant
xp = [SOC0_plant; ir0_plant];

% Pre-allocating storage
voltage_store_ff = zeros(Nsim,1);
u_store_ff = zeros(Nsim,1);
voltage_store_std = zeros(Nsim,1);
u_store_std = zeros(Nsim,1);

% Initializing 
xp_store_ff(:,1) = xp;
u_store_ff(1) = u0;
xp_store_std(:,1) = xp;
u_store_std(1) = u0;
[mpcData] = initMPCmodel(xp,mpcData);
% Initial voltage given initial SOC
v0 = OCVfromSOCtemp(xp(1),25,model) + mpcData.matrices.Cv*xp + mpcData.matrices.Dv*u0; 
xp = mpcData.matrices.A*xp + mpcData.matrices.B*u0;
xp_store_ff(:,2) = xp;
xp_store_std(:,2) = xp;
voltage_store_ff(1) = v0;
voltage_store_std(1) = v0;
mpcData_sf = mpcData;
mpcData_std = mpcData;
xp_ff = xp;  xp_std = xp;

for k = 1:Nsim
    
    % MPC     
    [uk_ff, mpcData_sf] = iterMPC(SOC_ref/100,xp_ff,mpcData_sf);  % Split future MPC
    [uk_std, mpcData_std] = iterMPC_standard(SOC_ref/100,xp_std,mpcData_std); % Standard MPC
    % Plant
    [voltage_ff, xp_ff] = iterModel(xp_ff,uk_ff,model,deltaT);
    [voltage_std, xp_std] = iterModel(xp_std,uk_std,model,deltaT);
    
    
    u_store_ff(k+1) = uk_ff; u_store_std(k+1) = uk_std;
    xp_store_ff(:,k+2) = xp_ff;   xp_store_std(:,k+2) = xp_std;
    voltage_store_ff(k+1) = voltage_ff;   voltage_store_std(k+1) = voltage_std;

end

toc


% time=time/60;     % Time in minutes
N = ones(1,length(time));

figure();
% plot current
subplot(311)
plot(0:time(end)-1,u_store_ff,'lineWidth',2); hold on;
plot(0:time(end)-1,u_store_std,'lineWidth',2);hold on;
plot(0:time(end)-1,mpcData.const.u_min*N,'--r','lineWidth',2);
xlim([0 (time(end)-1)])
grid on; title('Charging Current');xlabel('Time [sec]');ylabel('Current [A]');
legend('FF MPC','Standard MPC','Constraint');

% plot voltage
subplot(312)
plot(0:time(end)-1,voltage_store_ff,'lineWidth',2);hold on;
plot(0:time(end)-1,voltage_store_std,'lineWidth',2);hold on;
plot(0:time(end)-1,mpcData.const.v_max*N,'--r','lineWidth',2);
xlim([0 (time(end)-1)]);
ylim([voltage_store_ff(1) (mpcData.const.v_max+0.2)]); 
grid on; title('Voltage');xlabel('Time [sec]');ylabel('Voltage [V]');
legend('FF MPC','Standard MPC','Constraint');

% plot SOC
subplot(313)
plot(0:time(end)-1,xp_store_ff(1,1:end-1),'lineWidth',2);hold on;
plot(0:time(end)-1,xp_store_std(1,1:end-1),'lineWidth',2);hold on;
plot(0:time(end)-1,mpcData.const.z_max*N,'--r','lineWidth',2);
ylabel('SOC [%]');
grid on; title('State of Charge');xlabel('Time [sec]');
xlim([0 (time(end)-1)]);
legend('FF MPC','Standard MPC','Constraint');


