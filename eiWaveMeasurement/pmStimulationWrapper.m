

% to do --
% - get average of many runs
% - do overlapping sampling to increase sample rate of estimate

% Log Normal Inline
lnprm = @(mn,vr) [log((mn^2)/sqrt(vr+mn^2)), sqrt(log(vr/(mn^2)+1))];

% Time
tprm.dt = 0.01; % ms
tprm.T = 100; % ms


% Synaptic Parameters - alpha conductance
alpha = @(t, rise, fall) double(t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's...

exc.excRise = 0.3; % ms
exc.excFall = 5; % ms
exc.excAmp = 3e-9; % S
exc.excRev = 0e-3; % V
exc.excNumberMean = 30;  % Mean Value
exc.excNumberVar = 15; % Variance
exc.enPrms = lnprm(excNumberMean,excNumberVar); % Excitatory Number - LN Parameters
exc.excHardDelay = 5;
exc.excDelayMean = 2;
exc.excDelayVar = 3; 
exc.edPrms = lnprm(excDelayMean,excDelayVar); % Excitatory Delay - LN Parameters

inh.inhRise = 2.5; % ms
inh.inhFall = 10; % ms
inh.inhAmp = 3e-9; % S
inh.inhRev = -70e-3; % V
inh.inhDelay = 7; % ms
inh.inhNumberMean = 30; 
inh.inhNumberVar = 15; 
inh.inPrms = lnprm(inhNumberMean,inhNumberVar); % Inhibitory Number - LN Parameters
inh.inhHardDelay = 7;
inh.inhDelayMean = 2;
inh.inhDelayVar = 3;
inh.idPrms = lnprm(excDelayMean,excDelayVar); % Inhibitory Delay - LN Parameters

% Stimulation / Analysis Parameters
stim.vHold = -35e-3; % V
stim.modulationShift = 1; % 1/0 - do you randomize the phase?
stim.modulationDepth = 10e-3; %V
stim.modulationPeriod = 2; % ms
stim.noiseAmplitude = 15e-12; % A
stim.aCycles = 0.5; 
stim.method = 'simple'; % or "interleaved" or a numeric value to downsample interleaved

% Run Simulations
NR = 50;
results = cell(NR,1);
for nr = 1:NR
    if rem(nr*10,NR)==0
        msg = sprintf('%d/%d -- %d%%%%\\n',nr,NR,round(100*nr/NR));
        fprintf(1,msg);
    end
    results{nr} = pmStimulation(tprm,exc,inh,stim);
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


f=figure(1);
clf;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.7]);

subplot(3,1,1);
plot(base.tvec, 1e3*mean(holdVoltage,1),'k');
xlim([0 50]);
ylim([-45 -25]);
ylabel('mV');
title('hold voltage');
set(gca,'fontsize',16);

subplot(3,1,2);
hold on;
plot(base.tvec,mean(1e9*eConductance,1),'k');
plot(base.tvec,mean(1e9*iConductance,1),'r');
shadedErrorBar(base.tvec,mean(1e9*eConductance,1),std(1e9*eConductance,[],1),'k',0.5);
shadedErrorBar(base.tvec,mean(1e9*iConductance,1),std(1e9*iConductance,[],1),'r',0.5);
xlim([0 50]);
ylabel('nS');
legend('excitatory','inhibitory','location','northeast');
title('conductance distribution');
set(gca,'fontsize',16);

subplot(3,1,3);
hold on;
shadedErrorBar(base.tvec,mean(1e12*totalCurrent,1),std(1e12*totalCurrent,[],1),{'color','k','linewidth',1.5});
xlim([0 50]);
ylabel('pA');
title('Total Current Distribution');
set(gca,'fontsize',16);

%print(f,'-painters',fullfile(spath,'simulationPlots'),'-djpeg');

g = figure(2);
clf;
set(gcf,'units','normalized','outerposition',[0.55 0.1 0.4 0.7]);

subplot(3,2,1);
hold on;
plot(base.aWindowCenter, mean(1e9*cycConductance(1,:,:),3), 'color','k','linewidth',1.5);
plot(base.aWindowCenter, mean(1e9*estConductance(1,:,:),3), 'color','k','linestyle','none','linewidth',1.5,'marker','*');
xlim([0 50]);
ylim([-30 60]);
ylabel('nS');
legend('AvgConductance','Estimate');
title('Excitatory Conductance');
set(gca,'fontsize',16);

subplot(3,2,2);
hold on;
plot(base.aWindowCenter, mean(1e9*cycConductance(2,:,:),3), 'color','r','linewidth',1.5);
plot(base.aWindowCenter, mean(1e9*estConductance(2,:,:),3), 'color','r','linestyle','none','linewidth',1.5,'marker','*');
xlim([0 50]);
ylim([-30 60]);
ylabel('nS');
legend('AvgConductance','Estimate');
title('Inhibitory Conductance');
set(gca,'fontsize',16);

subplot(3,2,[3 4]);
hold on;
plot(base.aWindowCenter, mean(excDifference,2),'k');
plot(base.aWindowCenter, mean(inhDifference,2),'r');
shadedErrorBar(base.aWindowCenter, mean(excDifference,2), std(excDifference,[],2), 'k',0.5);
shadedErrorBar(base.aWindowCenter, mean(inhDifference,2), std(inhDifference,[],2), 'r',0.5);
xlim([0 50]);
ylim([-30 60]);
ylabel('nS');
title('Error in Conductance Estimation');
legend('Excitatory Error','Inhibitory Error');
set(gca,'fontsize',16);

subplot(3,2,5);
he = histogram(excr2,0:0.05:1);
he.FaceColor = 'k';
he.EdgeColor = 'k';
xlim([0 1]);
xlabel('R^{2}');
ylabel('counts');
title('Excitatory Estimation Accuracy');
set(gca,'fontsize',16);

subplot(3,2,6);
hi = histogram(inhr2,0:0.05:1);
hi.FaceColor = 'r';
hi.EdgeColor = 'r';
xlim([0 1]);
xlabel('R^{2}');
ylabel('counts');
title('Inhibitory Estimation Accuracy');
set(gca,'fontsize',16);

print(g,'-painters',fullfile(spath,'errorAnalysis-50x-PhaseShift'),'-djpeg');










    







