
% Time
tprm.dt = 0.01e-3; % s
tprm.T = 100e-3; % s

% Synaptic Parameters - alpha conductance 
alpha = @(t, rise, fall) double(t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's... 

exc.rise = 0.5e-3; % ms
exc.fall = 6e-3; % ms
exc.amp = 50e-12; % S
exc.rev = 0e-3; % V
exc.numberMean = 30;  % Mean Value 
exc.numberVar = 0; % Variance
exc.hardDelay = 0e-3; 
exc.delayMean = 5e-3; 
exc.delayVar = 0; 

inh.rise = 2e-3; % ms
inh.fall = 7e-3; % ms
inh.amp = 50e-12; % S
inh.rev = -70e-3; % V
inh.numberMean = 40; 
inh.numberVar = 0;
inh.hardDelay = 0e-3;
inh.delayMean = 5e-3;
inh.delayVar = 0e-3;

% Stimulation / Analysis Parameters 
stim.vHold = -35e-3; % V
stim.modulationShift = -1; % 1/0 - do you flip the phase? if -1, then randomize phase 
stim.modulationDepth = 15e-3; %V
stim.modulationPeriod = 2e-3; % s
stim.noiseAmplitude = 3e-12; % A
stim.aCycles = 1; 
stim.method = 'simple'; % or "interleaved" or a numeric value to downsample interleaved 

% Cell Parameters
cellprm.rs = 10e6; % access resistance (Ohms)
cellprm.rm = 150e6; % input resistance (Ohms)
cellprm.cm = 100e-12; % cell capacitance (F)
cellprm.em = -70e-3; % rest potential (V)

% Run simulation without open circuit
cellprm.openCircuit = 0; 
rtrue = vcStimulation(tprm,exc,inh,stim,cellprm);
% Run simulation with open circuit
cellprm.openCircuit = 1;
ropen = vcStimulation(tprm,exc,inh,stim,cellprm);
fprintf(1,'Simulations finished.\n');

% 







