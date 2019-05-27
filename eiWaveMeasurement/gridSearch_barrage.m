
% hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/eiWaveMeasurement';

% Time
dt = 0.05; % ms
T = 200; % ms
tvec = 0:dt:T; %vector (ms)
NT = length(tvec);

% Synaptic Parameters - alpha conductance
alpha = @(t, rise, fall) double(t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's...

excRise = 0.3; % ms
excFall = 5; % ms
excAmp = 3e-9; % S
excRev = 0e-3; % V
excFreq = 100; % 1/sec
excNumber = excFreq * T/1000; % number of inputs

inhRise = 2.5; % ms
inhFall = 10; % ms
inhAmp = 3e-9; % S
inhRev = -70e-3; % mV
inhFreq = 100; % 1/sec
inhNumber = inhFreq * T/1000; % Number of inputs

% Generate Conductances
eTrain = randi(T,excNumber, 1);
iTrain = randi(T,inhNumber, 1);
excConds = zeros(excNumber,NT);
inhConds = zeros(inhNumber,NT);
for ne = 1:excNumber
    cConductance = excAmp * alpha(tvec-eTrain(ne),excRise,excFall);
    cConductance(isnan(cConductance))=0;
    excConds(ne,:) = cConductance;
end
for ni = 1:inhNumber
    cConductance = inhAmp * alpha(tvec-iTrain(ni),inhRise,inhFall);
    cConductance(isnan(cConductance)) = 0;
    inhConds(ni,:) = cConductance;
end
eConductance = sum(excConds,1);
iConductance = sum(inhConds,1);

% Recording Parameters
capacitance = 30e-12; % Farads
inputResistance = 300e6; % Ohms
restPotential = -70e-3; % Volt

% Stimulation Parameters
modulationDepth = 10e-3; % mV
modulationPeriod = 0.5; % ms
holdVoltage = -35e-3 + modulationDepth/2*sin(2*pi*tvec/modulationPeriod); % Hold Voltage Command
eCurrent = eConductance .* (holdVoltage - excRev); % Excitatory Current
iCurrent = iConductance .* (holdVoltage - inhRev); % Inhibitory Current
synCurrent = eCurrent + iCurrent; % Total Synaptic Current
nCurrent = sqrt(15)*1e-12*randn(1,length(tvec)); % Noise current
totalCurrent = synCurrent + nCurrent; % Total Measured Current

% Analysis Cycles
aCycles = 2; % How many cycles of the modulation period to average over?
aWindowTime = aCycles * modulationPeriod; % Time of analysis window
aWindowSamples = aWindowTime/dt; % Samples of analysis window
aWindowStart = 1:aWindowSamples:NT; % Start of each window in samples
aWindowCenter = mean([tvec(aWindowStart(1:end-1));tvec(aWindowStart(2:end))],1); % Center of each window in time
NAW = length(aWindowStart); % Number of analysis windows
aResult = zeros(aWindowSamples,NAW-1,2); 
aLine = zeros(2,NAW-1);
estConductance = zeros(2,NAW-1); % Estimate of conductance [gExc; gInh]
residual = zeros(2,NAW-1); % Residual (based on average of each analysis window)
for naw = 1:NAW-1
    cSamples = aWindowStart(naw):aWindowStart(naw+1)-1;
    
    % Results --
    aResult(:,naw,1) = totalCurrent(cSamples);
    aResult(:,naw,2) = holdVoltage(cSamples);
    aLine(:,naw) = [aResult(:,naw,2) ones(aWindowSamples,1)] \ aResult(:,naw,1);
    
    % Convert to Conductances
    cgSyn = aLine(1,naw);
    ceSyn = -aLine(2,naw) / cgSyn;
    cgi = (cgSyn*ceSyn - cgSyn*excRev) / (inhRev - excRev);
    cge = cgSyn - cgi;
    estConductance(1,naw) = cge;
    estConductance(2,naw) = cgi;
    
    % Compute Residual
    cExcAverage = mean(eConductance(cSamples));
    cInhAverage = mean(iConductance(cSamples));
    residual(:,naw) = estConductance(:,naw) - [cExcAverage;cInhAverage];
end


% Plot Results
f = figure(1);
clf;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.7 0.8]);

subplot(2,2,1);
plot(tvec, 1e12*totalCurrent,'k','linewidth',1.5);
xlabel('Time (ms)');
ylabel('pA');
title('V-Clamp Recording');
set(gca,'fontsize',16);

subplot(2,2,3);
hold on;
plot(tvec, 1e12*eCurrent, 'k', 'linewidth',1.5);
plot(tvec, 1e12*iCurrent, 'r', 'linewidth',1.5);
legend('I_e','I_i','location','northeast');
xlabel('Time (ms)');
ylabel('pA');
title('Synaptic Currents');
set(gca,'fontsize',16);

subplot(2,2,2);
hold on;
% True Conductances
plot(tvec, 1e9*eConductance, 'k','linewidth',1.5);
plot(tvec, 1e9*iConductance, 'r','linewidth',1.5);
% Estimated Conductances
plot(aWindowCenter,1e9*estConductance(1,:),...
    'color','k','linestyle','none','marker','o','markerfacecolor','k','markersize',3);
plot(aWindowCenter,1e9*estConductance(2,:),...
    'color','r','linestyle','none','marker','o','markerfacecolor','r','markersize',3);
legend('G_e','G_i','location','northeast');
xlabel('Time (ms)');
ylabel('nS');
title('Synaptic Conductance w/ Estimates');
set(gca,'fontsize',16);

subplot(2,2,4);
hold on;
plot(aWindowCenter,1e9*residual(1,:),'color','k','linewidth',1.5);
plot(aWindowCenter,1e9*residual(2,:),'color','r','linewidth',1.5);
xlabel('Time (ms)');
ylabel('nS');
title('Residual Error');
set(gca,'fontsize',16);





















    


