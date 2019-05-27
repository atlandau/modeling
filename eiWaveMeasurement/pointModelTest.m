
hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/eiWaveMeasurement';

% Time
dt = 0.01; % ms
T = 40; % ms
tvec = 0:dt:T; %vector (ms)

% Synaptic Parameters - alpha conductance
alpha = @(t, rise, fall) (t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's...

excRise = 0.3; % ms
excFall = 5; % ms
excAmp = 2e-9; % S
excRev = 0e-3; % V
excDelay = 5; % ms

inhRise = 2.5; % ms
inhFall = 10; % ms
inhAmp = 2e-9; % S
inhRev = -70e-3; % mV
inhDelay = 5; % ms

% Total Conductance
eConductance = excAmp * alpha(tvec-excDelay, excRise, excFall);
iConductance = inhAmp * alpha(tvec-inhDelay, inhRise, inhFall);

% Recording Parameters
capacitance = 30e-12; % Farads
inputResistance = 300e6; % Ohms
restPotential = -70e-3; % Volt

% Stimulation Parameters
modulationDepth = 25e-3; % mV
modulationPeriod = 0.5;
holdVoltage = -35e-3 + modulationDepth/2*sin(2*pi*tvec/modulationPeriod);
eCurrent = eConductance .* (holdVoltage - excRev);
iCurrent = iConductance .* (holdVoltage - inhRev);









% Plot Results
f = figure(1);
clf;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.8]);

subplot(4,1,1);
plot(tvec, 1e3*holdVoltage, 'k','linewidth',1.5);
ylabel('V_m');
title('Holding Voltage');

subplot(4,1,2);
plot(tvec, 1e12*(eCurrent+iCurrent),'k','linewidth',1.5);
ylabel('pA');
title('V-Clamp Recording');

subplot(4,1,3);
hold on;
plot(tvec, 1e12*eCurrent, 'k', 'linewidth',1.5);
plot(tvec, 1e12*iCurrent, 'r', 'linewidth',1.5);
legend('I_e','I_i','location','northeast');
ylabel('pA');

subplot(4,1,4);
hold on;
plot(tvec, 1e9*eConductance, 'k','linewidth',1.5);
plot(tvec, 1e9*iConductance, 'r','linewidth',1.5);
legend('G_e','G_i','location','northeast');
ylabel('nS');


