

% to do --
% - get average of many runs
% - do overlapping sampling to increase sample rate of estimate

% Log Normal Inline
lnprm = @(mn,vr) [log((mn^2)/sqrt(vr+mn^2)), sqrt(log(vr/(mn^2)+1))];

% Time
tprm.dt = 0.01e-3; % s
tprm.T = 25e-3; % s


% Synaptic Parameters - alpha conductance
alpha = @(t, rise, fall) double(t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's...

exc.excRise = 0.3e-3; % ms
exc.excFall = 5e-3; % ms
exc.excAmp = 100e-12; % S
exc.excRev = 0e-3; % V
exc.excNumberMean = 50;  % Mean Value
exc.excNumberVar = 15; % Variance
exc.excHardDelay = 5e-3;
exc.excDelayMean = 2e-3;
exc.excDelayVar = 3e-3; 

inh.inhRise = 2.5e-3; % ms
inh.inhFall = 10e-3; % ms
inh.inhAmp = 150e-12; % S
inh.inhRev = -70e-3; % V
inh.inhNumberMean = 25; 
inh.inhNumberVar = 15;
inh.inhHardDelay = 7e-3;
inh.inhDelayMean = 2e-3;
inh.inhDelayVar = 3e-3;

% Stimulation / Analysis Parameters
stim.vHold = -35e-3; % V
stim.modulationShift = -1; % 1/0 - do you flip the phase?
stim.modulationDepth = 15e-3; %V
stim.modulationPeriod = 2e-3; % s
stim.noiseAmplitude = 15e-12; % A
stim.aCycles = 0.5; 
stim.method = 'simple'; % or "interleaved" or a numeric value to downsample interleaved

% Cell Parameters
cellprm.rs = 10e6; % access resistance (Ohms)
cellprm.rm = 500e6; % input resistance (Ohms)
cellprm.cm = 30e-12; % cell capacitance (F)
cellprm.em = -70e-3; % rest potential (V)

% Run Simulations
NR = 500;
results = cell(NR,1);
for nr = 1:NR
    if rem(nr*10,NR)==0
        msg = sprintf('%d/%d -- %d%%%%\\n',nr,NR,round(100*nr/NR));
        fprintf(1,msg);
    end
    results{nr} = vcStimulation(tprm,exc,inh,stim,cellprm);
end

%save(fullfile(spath,'results_noPhaseShiftMod.mat'),'results');


%% Plot Results
spath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/eiWaveMeasurement';

base = results{1};
holdVoltage = cell2mat(cellfun(@(c) c.holdVoltage, results, 'uni', 0));
eConductance = cell2mat(cellfun(@(c) c.eConductance, results, 'uni', 0));
iConductance = cell2mat(cellfun(@(c) c.iConductance, results, 'uni', 0));
totalCurrent = cell2mat(cellfun(@(c) c.totalCurrent, results, 'uni', 0));
estConductance = cell2mat(permute(cellfun(@(c) c.estConductance, results, 'uni', 0),[3 2 1]));
cycConductance = cell2mat(permute(cellfun(@(c) c.cycleConductance, results, 'uni', 0),[3 2 1]));
excDifference = permute(1e9*estConductance(1,:,:) - 1e9*cycConductance(1,:,:),[2 3 1]);
inhDifference = permute(1e9*estConductance(2,:,:) - 1e9*cycConductance(2,:,:),[2 3 1]);
excr2 = cellfun(@(c) c.excr2, results, 'uni', 1);
inhr2 = cellfun(@(c) c.inhr2, results, 'uni', 1);

xLim = [0 25];

% f=figure(1);
% clf;
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.7]);
% 
% subplot(3,1,1);
% plot(1e3*base.tvec, 1e3*mean(holdVoltage,1),'k');
% xlim(xLim);
% ylim([-45 -25]);
% ylabel('mV');
% title('hold voltage');
% set(gca,'fontsize',16);
% 
% subplot(3,1,2);
% hold on;
% plot(1e3*base.tvec,mean(1e9*eConductance,1),'k');
% plot(1e3*base.tvec,mean(1e9*iConductance,1),'r');
% shadedErrorBar(1e3*base.tvec,mean(1e9*eConductance,1),std(1e9*eConductance,[],1),'k',0.5);
% shadedErrorBar(1e3*base.tvec,mean(1e9*iConductance,1),std(1e9*iConductance,[],1),'r',0.5);
% xlim(xLim);
% ylabel('nS');
% legend('excitatory','inhibitory','location','northeast');
% title('conductance distribution');
% set(gca,'fontsize',16);
% 
% subplot(3,1,3);
% hold on;
% shadedErrorBar(1e3*base.tvec,mean(1e12*totalCurrent,1),std(1e12*totalCurrent,[],1),{'color','k','linewidth',1.5});
% xlim(xLim);
% ylabel('pA');
% title('Total Current Distribution');
% set(gca,'fontsize',16);

%print(f,'-painters',fullfile(spath,'simulationPlots'),'-djpeg');

g = callFigs(2);
clf;
% set(gcf,'units','normalized','outerposition',[0.55 0.1 0.4 0.7]);

subplot(3,2,1);
hold on;
plot(1e3*base.aWindowCenter, mean(1e9*cycConductance(1,:,:),3), 'color','k','linewidth',1.5);
plot(1e3*base.aWindowCenter, mean(1e9*estConductance(1,:,:),3), 'color','k','linestyle','none','linewidth',1.5,'marker','*');
xlim(xLim);
ylim([-1 3.5]);
ylabel('nS');
legend('AvgConductance','Estimate');
title('Excitatory Conductance');
set(gca,'fontsize',16);

subplot(3,2,2);
hold on;
plot(1e3*base.aWindowCenter, mean(1e9*cycConductance(2,:,:),3), 'color','r','linewidth',1.5);
plot(1e3*base.aWindowCenter, mean(1e9*estConductance(2,:,:),3), 'color','r','linestyle','none','linewidth',1.5,'marker','*');
xlim(xLim);
ylim([-1 3.5]);
ylabel('nS');
legend('AvgConductance','Estimate');
title('Inhibitory Conductance');
set(gca,'fontsize',16);

subplot(3,2,[3 4]);
hold on;
plot(1e3*base.aWindowCenter, mean(excDifference,2),'k');
plot(1e3*base.aWindowCenter, mean(inhDifference,2),'r');
shadedErrorBar(1e3*base.aWindowCenter, mean(excDifference,2), std(excDifference,[],2), 'k',0.5);
shadedErrorBar(1e3*base.aWindowCenter, mean(inhDifference,2), std(inhDifference,[],2), 'r',0.5);
xlim(xLim);
ylim([-3.5 3.5]);
ylabel('nS');
title('Error in Conductance Estimation');
legend('Excitatory Error','Inhibitory Error');
set(gca,'fontsize',16);

subplot(3,2,5);
he = histogram(excr2,0:0.025:1);
he.FaceColor = 'k';
he.EdgeColor = 'k';
xlim([0 1]);
xlabel('R^{2}');
ylabel('counts');
title('Excitatory Estimation Accuracy');
set(gca,'fontsize',16);

subplot(3,2,6);
hi = histogram(inhr2,0:0.025:1);
hi.FaceColor = 'r';
hi.EdgeColor = 'r';
xlim([0 1]);
xlabel('R^{2}');
ylabel('counts');
title('Inhibitory Estimation Accuracy');
set(gca,'fontsize',16);

% print(g,'-painters',fullfile(spath,'errorAnalysis-50x-PhaseShift'),'-djpeg');










    







