function o2_dynamicKappa_bufferDynamics011(konIDX,kdIDX,ampIDX)

% Timing Parameters
dt = 1e-6; % seconds
ds = 1000; % downsample factor
T = 0.1;

% General Biophysical Parameters
ks = 20; % endogenous binding ratio -- assuming instant in this model
p.beta = ks*1000/12; % modeled from Sabatini et al., 2002
p.rest = 75e-9; % Resting calcium (M)
p.v = (4/3) * pi * 0.4^3 * 1e-15; % Volume spine (L)
p.bconc = 300e-6; % Concentration exogenous buffer (M)

% Grid parameters
kon = [1e6 5e6 1e7 5e7 1e8 5e8 1e9 5e9 1e10]; % Kon values
kd = [400e-9 800e-9 2e-6 5e-6 10e-6 40e-6 100e-6]; % Kds
amp = [0.5e-12 1e-12 2e-12 5e-12 10e-12 20e-12 50e-12]; % Amplitude values

p.kon = kon(konIDX); % Kon
p.koff = kon(konIDX) * kd(kdIDX); % Koff
p.amplitude = amp(ampIDX); % Peak amplitude calcium current

% Useful inlines
kappa = @(p,ca) p.bconc / (ca + p.koff/p.kon); % Inline for equilibrium binding ratio
ica = @(t,amp) (1/ks) * amp * exp(-(t-2e-3).^2./(0.55e-3).^2) .* (t >= 1e-3) .* (t <= 3e-3); 

currentKappa = kappa(p,p.rest);
iState = [p.rest p.rest*currentKappa];

y = eulerapp(@(t,y) bufferDynamics_011(t,y,ica,p),[0 T],iState,dt,ds); % numerical approximation of calcium and buffer

% Make grid parameter structure so we know what the values were
gridParameters.kon = kon;
gridParameters.kd = kd;
gridParameters.amp = amp;

% Save y (data), p (prms), gridParameters (current search)
name = sprintf('results_Kon%d_Kd%d_Amp%d.mat',konIDX,kdIDX,ampIDX);
savepath = fullfile('~/bufferCapacity_190523',name);
save(savepath,'y','p','gridParameters');





