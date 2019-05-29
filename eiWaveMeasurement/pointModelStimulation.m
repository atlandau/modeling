


% to do --
% - get average of many runs
% - do overlapping sampling to increase sample rate of estimate

% Log Normal Inline
lnprm = @(mn,vr) [log((mn^2)/sqrt(vr+mn^2)), sqrt(log(vr/(mn^2)+1))];

% Time
dt = 0.01; % ms
T = 100; % ms
tvec = 0:dt:T; %vector (ms)
NT = length(tvec);
stimTime = 5; % ms
stimSample = find(tvec>=stimTime,1); % sample of stimulation time

% Synaptic Parameters - alpha conductance
alpha = @(t, rise, fall) double(t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's...

excRise = 0.3; % ms
excFall = 5; % ms
excAmp = 3e-9; % S
excRev = 0e-3; % V
excNumberMean = 30;  % Mean Value
excNumberVar = 15; % Variance
enPrms = lnprm(excNumberMean,excNumberVar); % Excitatory Number - LN Parameters
excHardDelay = 5;
excDelayMean = 2;
excDelayVar = 3; 
edPrms = lnprm(excDelayMean,excDelayVar); % Excitatory Delay - LN Parameters

inhRise = 2.5; % ms
inhFall = 10; % ms
inhAmp = 3e-9; % S
inhRev = -70e-3; % mV
inhDelay = 7; % ms
inhNumberMean = 30; 
inhNumberVar = 15; 
inPrms = lnprm(inhNumberMean,inhNumberVar); % Inhibitory Number - LN Parameters
inhHardDelay = 7;
inhDelayMean = 2;
inhDelayVar = 3;
idPrms = lnprm(excDelayMean,excDelayVar); % Inhibitory Delay - LN Parameters


% Generate Conductances
eNumber = round(lognrnd(enPrms(1),enPrms(2)));
iNumber = round(lognrnd(inPrms(1),inPrms(2)));
% Generate timing trains, resolve to dt, add stimulation time as hard delay
eTrain = excHardDelay + round(lognrnd(edPrms(1),edPrms(2),[eNumber 1])/dt)*dt;
iTrain = inhHardDelay + round(lognrnd(idPrms(1),idPrms(2),[iNumber 1])/dt)*dt;
excConds = zeros(eNumber,NT);
inhConds = zeros(iNumber,NT);
for ne = 1:eNumber
    cConductance = excAmp * alpha(tvec-eTrain(ne),excRise,excFall);
    cConductance(isnan(cConductance))=0;
    excConds(ne,:) = cConductance;
end
for ni = 1:iNumber
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
modulationPeriod = 2; % ms
holdVoltage = -35e-3 + modulationDepth/2*sin(2*pi*tvec/modulationPeriod); % Hold Voltage Command
eCurrent = eConductance .* (holdVoltage - excRev); % Excitatory Current
iCurrent = iConductance .* (holdVoltage - inhRev); % Inhibitory Current
synCurrent = eCurrent + iCurrent; % Total Synaptic Current
nCurrent = sqrt(15)*1e-12*randn(1,length(tvec)); % Noise current
totalCurrent = synCurrent + nCurrent; % Total Measured Current

% Analysis Cycles
aCycles = 0.5; % How many cycles of the modulation period to average over?
aWindowTime = aCycles * modulationPeriod; % Time of analysis window
aWindowSamples = aWindowTime/dt; % Samples of analysis window
aWindowStart = 1:aWindowSamples:NT; % Start of each window in samples
aWindowCenter = mean([tvec(aWindowStart(1:end-1));tvec(aWindowStart(2:end))],1); % Center of each window in time
NAW = length(aWindowStart); % Number of analysis windows
aResult = zeros(aWindowSamples,NAW-1,2); 
aLine = zeros(2,NAW-1);
estConductance = zeros(2,NAW-1); % Estimate of conductance [gExc; gInh]
cycleConductance = zeros(2,NAW-1); % Average of conductance per analysis cycle
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
    cycleConductance(:,naw) = [mean(eConductance(cSamples)); mean(iConductance(cSamples))];
    residual(:,naw) = estConductance(:,naw) - cycleConductance(:,naw);
end

rmsError = sqrt(mean((1e9*residual).^2,2)); % rms error
[~,excGOF] = fit(1e9*cycleConductance(1,:)',1e9*estConductance(1,:)','poly1'); % linear fit to get R^2
[~,inhGOF] = fit(1e9*cycleConductance(2,:)',1e9*estConductance(2,:)','poly1'); % linear fit to get R^2
excr2 = excGOF.rsquare;
inhr2 = inhGOF.rsquare;




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
bar([rmsError' excr2 inhr2],'FaceColor','k');
line(xlim,[1 1],'color','k');
set(gca,'xtick',1:4);
set(gca,'xticklabel',{'RMS(Res)-Exc','RMS(Res)-Inh','Exc R^{2}','Inh R^{2}'});
set(gca,'xticklabelrotation',45);
ylabel('nS      or      R^{2}');
title('Summary Stats');
set(gca,'fontsize',16);


















    







