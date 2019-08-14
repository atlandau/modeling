
% Time
tprm.dt = 0.01e-3; % s
tprm.T = 40e-3; % s

% Synaptic Parameters - alpha conductance 
alpha = @(t, rise, fall) double(t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's... 

exc.rise = 0.5e-3; % ms
exc.fall = 6e-3; % ms
exc.amp = 50e-12; % S
exc.rev = 0e-3; % V
exc.numberMean = 80;  % Mean Value 
exc.numberVar = 0; % Variance
exc.hardDelay = 0e-3; 
exc.delayMean = 1e-3; 
exc.delayVar = 0; 

inh.rise = 2e-3; % ms
inh.fall = 7e-3; % ms
inh.amp = 50e-12; % S
inh.rev = -70e-3; % V
inh.numberMean = 100; 
inh.numberVar = 0;
inh.hardDelay = 0e-3;
inh.delayMean = 3e-3;
inh.delayVar = 0e-3;

% Stimulation / Analysis Parameters 
stim.vHold = -35e-3; % V
stim.modulationShift = 0; % 1/0 - do you flip the phase? if -1, then randomize phase 
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
result = vcStimulation(tprm,exc,inh,stim,cellprm);
% Run simulation with open circuit
cellprm.openCircuit = 1;
ropen = vcStimulation(tprm,exc,inh,stim,cellprm);
fprintf(1,'Simulations finished.\n');

%% Plots

tvec = 1e3*(result.tvec - tprm.T/2);
xLim = [-2 1e3*tprm.T/2];

NAW = size(result.aLine,2);
awidx = floor(NAW/2)+1:NAW;
awTime = 1e3*(result.aWindowCenter(awidx)-tprm.T/2);

figure(1);
clf;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.7]);

subplot(2,2,1); hold on;
plot(awTime,1e9*result.cycleConductance(1,awidx),'k','linewidth',1.5);
plot(awTime,1e9*result.cycleConductance(2,awidx),'r','linewidth',1.5);
plot(awTime,1e9*result.estConductance(1,awidx),'k','marker','*','linestyle','none');
plot(awTime,1e9*result.estConductance(2,awidx),'r','marker','*','linestyle','none');
xlim(xLim);
ylabel('nS');
title('True Result');
set(gca,'fontsize',16);

subplot(2,2,3); hold on;
plot(awTime,1e9*ropen.cycleConductance(1,awidx),'k','linewidth',1.5);
plot(awTime,1e9*ropen.cycleConductance(2,awidx),'r','linewidth',1.5);
plot(awTime,1e9*ropen.estConductance(1,awidx),'k','marker','*','linestyle','none');
plot(awTime,1e9*ropen.estConductance(2,awidx),'r','marker','*','linestyle','none');
xlim(xLim);
ylabel('nS');
title('What-If Result');
set(gca,'fontsize',16);

vLim = [-40 -35];
subplot(2,2,2); hold on;
cmap = varycolor(NAW/2);
for naw = 1:length(awidx), cidx = awidx(naw);
    plot(1e3*result.aResult(:,cidx,2),1e12*result.aResult(:,cidx,1),'color',cmap(naw,:),'marker','*','linestyle','none');
    plot(1e3*ropen.aResult(:,cidx,2),1e12*ropen.aResult(:,cidx,1),'color',cmap(naw,:),'linewidth',1,'marker','*');
end
xlim(vLim);
xlabel('mV');
ylabel('pA');

subplot(2,2,4); hold on;
for naw = 1:length(awidx), cidx = awidx(naw);
    line(vLim, 1e9*(result.aLine(1,cidx)-ropen.aLine(1,cidx))*vLim + 1e12*(result.aLine(2,cidx)-ropen.aLine(2,cidx)),'color',cmap(naw,:));
%     
%     line(vLim, 1e9*result.aLine(1,cidx)*vLim + 1e12*result.aLine(2,cidx),'color',cmap(naw,:),'marker','*');
%     line(vLim, 1e9*ropen.aLine(1,cidx)*vLim + 1e12*ropen.aLine(2,cidx),'color',cmap(naw,:));
end
xlim(vLim);

    



