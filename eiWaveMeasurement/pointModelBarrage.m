
hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/eiWaveMeasurement';

% Time
dt = 0.05; % ms
T = 1000; % ms
tvec = 0:dt:T; %vector (ms)

% Synaptic Parameters - alpha conductance
alpha = @(t, rise, fall) (t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's...

excRise = 0.3; % ms
excFall = 5; % ms
excAmp = 3e-9; % S
excRev = 0e-3; % V
excFreq = 50; % 1/sec

inhRise = 2.5; % ms
inhFall = 10; % ms
inhAmp = 3e-9; % S
inhRev = -70e-3; % mV
inhFreq = 50; % 1/sec

% Total Conductance
eConductance = excAmp * alpha(tvec-excDelay, excRise, excFall);
iConductance = inhAmp * alpha(tvec-inhDelay, inhRise, inhFall);

% Recording Parameters
capacitance = 30e-12; % Farads
inputResistance = 300e6; % Ohms
restPotential = -70e-3; % Volt

% Stimulation Parameters
modulationDepth = 10e-3; % mV
modulationPeriod = 0.5;
holdVoltage = -35e-3 + modulationDepth/2*sin(2*pi*tvec/modulationPeriod);
eCurrent = eConductance .* (holdVoltage - excRev);
iCurrent = iConductance .* (holdVoltage - inhRev);
synCurrent = eCurrent + iCurrent;
nCurrent = sqrt(10)*1e-12*randn(1,length(tvec));
totalCurrent = synCurrent + nCurrent;

% Analysis Cycles
aCycles = 4;
aWindowTime = aCycles * modulationPeriod;
aWindowSamples = aWindowTime/dt;
aWindowStart = 1:aWindowSamples:length(tvec);
aWindowCenter = mean([tvec(aWindowStart(1:end-1));tvec(aWindowStart(2:end))],1);
NAW = length(aWindowStart);
aResult = zeros(aWindowSamples,NAW-1,2);
aLine = zeros(2,NAW-1);
estConductance = zeros(2,NAW-1);
for naw = 1:NAW-1
    % Results --
    aResult(:,naw,1) = totalCurrent(aWindowStart(naw):aWindowStart(naw+1)-1);
    aResult(:,naw,2) = holdVoltage(aWindowStart(naw):aWindowStart(naw+1)-1);
    aLine(:,naw) = [aResult(:,naw,2) ones(aWindowSamples,1)] \ aResult(:,naw,1);
    
    % Convert to Conductances
    cgSyn = aLine(1,naw);
    ceSyn = -aLine(2,naw) / cgSyn;
    cgi = (cgSyn*ceSyn - cgSyn*excRev) / (inhRev - excRev);
    cge = cgSyn - cgi;
    estConductance(1,naw) = cge;
    estConductance(2,naw) = cgi;
end


% Plot Results
f = figure(1);
clf;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.8]);

subplot(4,1,1);
plot(tvec, 1e12*totalCurrent,'k','linewidth',1.5);
ylabel('pA');
title('V-Clamp Recording');

subplot(4,1,2);
hold on;
plot(tvec, 1e12*eCurrent, 'k', 'linewidth',1.5);
plot(tvec, 1e12*iCurrent, 'r', 'linewidth',1.5);
legend('I_e','I_i','location','northeast');
ylabel('pA');

subplot(4,1,3);
hold on;
% True Conductances
plot(tvec, 1e9*eConductance, 'k','linewidth',1.5);
plot(tvec, 1e9*iConductance, 'r','linewidth',1.5);
% Estimated Conductances
plot(aWindowCenter,1e9*estConductance(1,:),'color','k','linestyle','none','marker','o','markerfacecolor','k');
plot(aWindowCenter,1e9*estConductance(2,:),'color','r','linestyle','none','marker','o','markerfacecolor','r');
legend('G_e','G_i','location','northeast');
ylabel('nS');


subplot(4,1,4);
hold on;
cmap = varycolor(NAW-1);
for naw = 1:NAW-1
    plot(aResult(:,naw,2),aResult(:,naw,1),...
        'linestyle','none','color',cmap(naw,:),'marker','o','markerfacecolor',cmap(naw,:));
    line([-80e-3 10e-3],aLine(1,naw)*[-80e-3 10e-3]+aLine(2,naw),'color',cmap(naw,:),'linewidth',1.5);
end
set(gca,'xticklabel',1e3*get(gca,'xtick'));
set(gca,'yticklabel',1e12*get(gca,'ytick'));
xlim(1e-3*[-80 10]);


























    


