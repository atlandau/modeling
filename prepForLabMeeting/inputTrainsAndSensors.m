
% - Fluorescent Buffer Parameteres come directly from Saba, 2002
% - This is where I get values for Kd, 1/amp, tau, and [FB]
% - Assuming all are diffusion limited (5e8/M/s)
% - *** will eventually include a Mg term *** 
% The fluo-5f tau and peak are guesses by measuring shit in fiji

hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/prepForLabMeeting';
lmPath = fullfile(cdSL,'presentations/LabMeeting/190906');
dacPath = '/Users/LandauLand/Documents/Research/SabatiniLab/presentations/DACs/DAC2';

twoPanelSize = [0.29 0.27 0.21 0.38];
threePanelSize = [0.29 0.27 0.21*2/3 0.38*2/3];


%% Choose optimal Bt, Kd, Kon, then look at traces evoked by AP trains with different frequency

% ODE Parameters
tspan = [0,0.2];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.001);
dt = 0.0001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 0.5e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 50e-9; % 50nM
amplitude = 2e-12; % 7pA
sysProperties = [spineVolume, extrusionRate, restConcentration, amplitude];

% Endogenous Buffer
endogenousKon = 5e8;
endogenousKd = 175e-6;
endogenousBt = endogenousKd*20;
bufferProperties.onRate = endogenousKon;
bufferProperties.kd = endogenousKd;
bufferProperties.totalConcentration = endogenousBt;
kappa = bufferProperties.totalConcentration / (restConcentration + bufferProperties.kd);

iState = restConcentration*[1 kappa];

% AP Properties
frequency = 1000./fliplr([5 7 10 15 25 50 75 100 125 150 175 200]);
NF = length(frequency);
numAPs = [1 2 3 4 5];
NAP = length(numAPs);

msg = '';
sols = cell(1,NF,NAP);
values = zeros(NT,2,NF,NAP);
for nap = 1:NAP
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf('working... numAPs %d/%d\n',nap,NAP);
    fprintf(1,msg);
    for nf = 1:NF
        onsetTimes = linspace(0, (numAPs(nap)-1)/frequency(nf), numAPs(nap));

        prmODE = @(t,y) massAction_wTrain(t,y,sysProperties,bufferProperties,onsetTimes);
        sols{1,nf,nap} = ode23s(prmODE,tspan,iState,odeOptions);
        values(:,:,nf,nap) = deval(sols{1,nf,nap},tvec)';
    end
end
pkCalcium = squeeze(max(values(:,1,:,:),[],1));
relCalcium = pkCalcium ./ min(pkCalcium(:));


figure(189); clf;
set(gcf,'units','normalized','outerposition',[0.29 0.27 0.21 0.38]);
hold on;
cmap = flipud([linspace(0.1,1,NF)'.^(1/2),zeros(NF,2)]);
for nf = NF:-1:1
    plot(1e3*tvec,1e6*values(:,1,nf,2),'color',cmap(nf,:),'linewidth',1.5);
end
xlabel('Time (ms)');
ylabel('[Ca^{2+}] (然)');
title('Frequency Change w/ 2 APs');
set(gca,'fontsize',18);

figure(190); clf;
set(gcf,'units','normalized','outerposition',[0.6 0.27 0.21 0.38]);
hold on;
cmap = [linspace(0.1,1,NF)',zeros(NF,2)];
for nf = 1:NF
    plot(numAPs,1e6*pkCalcium(nf,:),'color',cmap(nf,:),'linewidth',1.5,'marker','*');
end
xlim([0.5 NAP+0.5]);
% ylim([0 max(relCalcium(:))*1.1]);
set(gca,'xtick',1:NAP);
set(gca,'xticklabel',numAPs);
% set(gca,'ytick',1:max(relCalcium(:)));
cb = colorbar('Location','Manual','Position',[0.25 0.45 0.06 0.45]);
colormap(gca,cmap);
cb.Ticks = [0 1];
cb.TickLabels = 1000./frequency([1 end]);
ylabel(cb,'Period b/w Spikes (ms)');

xlabel('Number of APs');
ylabel('Peak Calcium Influx (然)');
set(gca,'fontsize',16);



%% How does calmodulin respond to a train above the critical frequency for accumulation?

% ODE Parameters
tspan = [0 2];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.001);
dt = 0.0001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

% System Properties
spineVolume = 0.5e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 50e-9; % 50nM
amplitude = 2e-12; % 7pA
sysProperties = [spineVolume, extrusionRate, restConcentration, amplitude];

% Endogenous Buffer
bufferProperties.onRate = 5e8;
bufferProperties.kd = 175e-6;
bufferProperties.totalConcentration = bufferProperties.kd*20;
kappa = bufferProperties.totalConcentration / (restConcentration + bufferProperties.kd);

% CaM Properties
camConcentration = 1e-6;
camProperties = getCaMProps();
iStateCaM = getRestingCaM(camConcentration,restConcentration);

iState = [restConcentration; restConcentration*kappa; iStateCaM];

% AP Properties
frequency = 1000./fliplr([3 5 7 10 15 25 50 100 125 200 10000]);
NF = length(frequency);

msg = '';
sols = cell(1,NF);
values = zeros(NT,11,NF);
for nf = 1:NF
    cPeriod = 1/frequency(nf);
    onsetTimes = 0:cPeriod:cPeriod*1;%tspan(2);

    prmODE = @(t,y) massActionCaM_wTrain(t,y,sysProperties,bufferProperties,camProperties,onsetTimes);
    sols{1,nf} = ode23s(prmODE,tspan,iState,odeOptions);
    values(:,:,nf) = deval(sols{1,nf},tvec)';
end
pkCalcium = squeeze(max(values(:,1,:),[],1));
relCalcium = pkCalcium ./ min(pkCalcium(:));

pkCaM22 = squeeze(max(values(:,end,:),[],1));

%%
figure(189); clf;
set(gcf,'units','normalized','outerposition',[0.05 0.27 0.21 0.38]);
hold on;
cmap = flipud([linspace(0.1,1,NF)'.^(1/2),zeros(NF,2)]);
for nf = NF:-1:1
    plot(1e3*tvec,1e6*values(:,1,nf),'color',cmap(nf,:),'linewidth',1.5);
end
xlim([0 80]);
xlabel('Time (ms)');
ylabel('[Ca^{2+}] (然)');
title('Calcium Influx Train APs');
set(gca,'fontsize',18);
% print(gcf,'-painters',fullfile(dacPath,'CaM-CalciumTraces'),'-djpeg');

figure(190); clf;
set(gcf,'units','normalized','outerposition',[0.26 0.27 0.21 0.38]);
hold on;
for nf = NF:-1:1
    plot(1e3*tvec,1e6*values(:,end,nf),'color',cmap(nf,:),'linewidth',1.5);
end
xlim([0 80]);
xlabel('Time (ms)');
ylabel('[CaM]_{N2C2} (然)');
title('CaM Train APs');
set(gca,'fontsize',18);
% print(gcf,'-painters',fullfile(dacPath,'CaM-Traces'),'-djpeg');

figure(191); clf;
set(gcf,'units','normalized','outerposition',[0.47 0.27 0.21 0.38]);
hold on;
plot(frequency, 100*pkCaM22/camConcentration, 'color','k','marker','o','linewidth',1.5);
xlabel('AP Frequency');
ylabel('%CaM N2-C2 Activation');
set(gca,'ytick',0:0.01:0.02);
set(gca,'fontsize',18);
inset = axes('Position',[0.5 0.2 0.35 0.35]);
box(inset,'on');
hold(inset,'on');
plot(pkCalcium/min(pkCalcium),pkCaM22/min(pkCaM22),'marker','o','color','k','markerfacecolor','k');
plot(0:10,0:10,'linestyle','--','color','k','linewidth',0.5);
xlim([0 10]);
ylim([0 10]);
set(gca,'xtick',0:3:9);
set(gca,'ytick',0:3:9);
xlabel('Peak Calcium');
ylabel('CaM N2-C2 Activation');
set(gca,'fontsize',14);
print(gcf,'-painters',fullfile(dacPath,'CaM-PercentActivation'),'-djpeg');


figure(192); clf;
set(gcf,'units','normalized','outerposition',[0.68 0.27 0.21 0.38]);
hold on;
plot(pkCalcium,pkCaM22,'linestyle','none','marker','o','color','k','markerfacecolor','k');
xlabel('Peak Calcium');
ylabel('CaM N2-C2 Activation');
set(gca,'fontsize',18);

%%
figure(190); clf;
set(gcf,'units','normalized','outerposition',[0.6 0.27 0.21 0.38]);
hold on;
cmap = [linspace(0.1,1,NF)',zeros(NF,2)];
for nf = 1:NF
    plot(numAPs,1e6*pkCalcium(nf,:),'color',cmap(nf,:),'linewidth',1.5,'marker','*');
end
xlim([0.5 NAP+0.5]);
% ylim([0 max(relCalcium(:))*1.1]);
set(gca,'xtick',1:NAP);
set(gca,'xticklabel',numAPs);
% set(gca,'ytick',1:max(relCalcium(:)));
cb = colorbar('Location','Manual','Position',[0.25 0.45 0.06 0.45]);
colormap(gca,cmap);
cb.Ticks = [0 1];
cb.TickLabels = 1000./frequency([1 end]);
ylabel(cb,'Period b/w Spikes (ms)');

xlabel('Number of APs');
ylabel('Peak Calcium Influx (然)');
set(gca,'fontsize',16);














%%

% use the best alpha fit from the glutamate binding profile to make an
% EPSP-AP pairing input and then put the AP at different times

% we can even use your measurements to amplify NMDAR influx accordingly

% 
























%% ODE for 1 spine with mass action buffer and AP train
function dydt = massAction_wTrain(t,y,sysProperties,bufferProperties,onsetTimes)
    freeCalcium = y(1); 
    boundBuffer = y(2); 
    freeBuffer = bufferProperties.totalConcentration - boundBuffer;
    
    % Calcium current stuff
    restCurrent = sysProperties(2) * sysProperties(3);
    curr = restCurrent + icaTrain(t,sysProperties(4),onsetTimes)/(2*96485*sysProperties(1));
    extrusion = sysProperties(2) * freeCalcium; 

    % Fluorescent Buffer Reaction
    fluorAssociation = freeCalcium * freeBuffer * bufferProperties.onRate; 
    fluorDissociation = boundBuffer * (bufferProperties.onRate * bufferProperties.kd);
    
    %Output
    dydt(1,1) = curr - extrusion - fluorAssociation + fluorDissociation;
    dydt(2,1) = fluorAssociation - fluorDissociation;
end

%% ODE for mass action buffer, CaM, and AP train
function dydt = massActionCaM_wTrain(t,y,sysProperties,bufferProperties,camProperties,onsetTimes)
    freeCalcium = y(1);
    boundBuffer = y(2);
    freeBuffer = bufferProperties.totalConcentration - boundBuffer;
    
    % Calcium current stuff
    restCurrent = sysProperties(2) * sysProperties(3);
    curr = restCurrent + icaTrain(t,sysProperties(4),onsetTimes)/(2*96485*sysProperties(1));
    extrusion = sysProperties(2) * freeCalcium; 

    % Fluorescent Buffer Reaction
    bufferAssociation = freeCalcium * freeBuffer * bufferProperties.onRate; 
    bufferDissociation = boundBuffer * (bufferProperties.onRate * bufferProperties.kd);
    
    % CaM On/Off Reactions
    C0_on = 2*camProperties.konTC * y(1) * y(3:5); %[00,10,20]
    C1_on = camProperties.konRC * y(1) * y(6:8); %[01,11,21]
    C1_off = camProperties.koffTC * y(6:8); %[01,11,21]
    C2_off = 2*camProperties.koffRC * y(9:11); %[02,12,22]
    N0_on = 2*camProperties.konTN * y(1) * y([3 6 9]); %[00,01,02]
    N1_on = camProperties.konRN * y(1) * y([4 7 10]); %[10,11,12]
    N1_off = camProperties.koffTN * y([4 7 10]); %[10,11,12]
    N2_off = 2*camProperties.koffRN * y([5 8 11]); %[20,21,22]
    
    % Transition Matrix
    cLobeTransitions = [-C0_on + C1_off, C0_on - C1_off - C1_on + C2_off, C1_on - C2_off];
    nLobeTransitions = [-N0_on + N1_off, N0_on - N1_off - N1_on + N2_off, N1_on - N2_off]';
    transitionMatrix = cLobeTransitions + nLobeTransitions;
    
    % Calcium Movement
    bufferExchange = bufferAssociation - bufferDissociation;
    camExchange = sum(C0_on + C1_on + N0_on + N1_on - C1_off - C2_off - N1_off - N2_off);
    
    %Output
    dydt(1,1) = curr - extrusion - bufferExchange - camExchange;
    dydt(2,1) = bufferExchange;
    dydt(3:11,1) = transitionMatrix(:);
end
































































