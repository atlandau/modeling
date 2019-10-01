
% - Fluorescent Buffer Parameteres come directly from Saba, 2002
% - This is where I get values for Kd, 1/amp, tau, and [FB]
% - Assuming all are diffusion limited (5e8/M/s)
% - *** will eventually include a Mg term *** 
% The fluo-5f tau and peak are guesses by measuring shit in fiji

hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/caBuffering_RecreateExperimentalData';
fpath = fullfile(hpath,'figures');

% Spine Parameters
systemPrms.beta=1400.0;
systemPrms.rest=5.0e-8;
systemPrms.ks=20.0;
systemPrms.v=(4/3)*pi*0.7^3*1e-15;
iCalciumAmp = 21e-12; 
systemPrms.amp = iCalciumAmp;

u0 = systemPrms.rest;
tspan = [0,0.5];
dt = 0.000001;
tvec = tspan(1):dt:tspan(2);

% Instant Solution
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',dt);
[trueTime,trueSol] = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms),tspan,u0,odeOptions);
plot(trueTime,trueSol);


%% Optimize endogenous:[Kon,Kd,Bt], rest, beta, and iCalcium to get following results with fluorophores in the cell

% ODE Parameters
tspan = [0,0.2];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.02);
dt = 0.00001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 1e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 70e-9; % 70nM
amplitude = 7e-12; % 7pA

% Properties of OGB
numOnRate = 75;
onRateRange = [4 12];
onRate = logspace(onRateRange(1),onRateRange(2),numOnRate); % 1/M/s
ogbConcentration = 50e-6;% 50然 example
ogbKD = 205e-9; % 205nM

fluorProperties.totalConcentration = ogbConcentration;
fluorProperties.kd = ogbKD;

fluorKappa = ogbConcentration / (restConcentration + ogbKD);
initialState = restConcentration * [1 fluorKappa];

sols = cell(1,numOnRate);
values = zeros(NT,2,numOnRate); 
for i = 1:numOnRate
    fluorProperties.onRate = onRate(i);
    
    prmODE = @(t,y) justFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,fluorProperties);
    sols{i} = ode23s(prmODE, tspan, initialState, odeOptions);
    
    values(:,:,i) = deval(sols{i},tvec)';
end

pkCalcium = max(squeeze(values(:,1,:)),[],1); % Max of calcium transient

plotIdx = 1:3:numOnRate;
rMap = [linspace(0.2,1,length(plotIdx))',zeros(length(plotIdx),2)];
figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.2 0.5 0.7 0.4]);
set(gcf,'DefaultAxesColorOrder',rMap);

subplot(1,3,1);
plot(1000*tvec, 1e6*squeeze(values(:,1,plotIdx)), 'linewidth',1.5);
xlim([0 25]);
yLimCalcium = ylim;
xlabel('Time (ms)');
ylabel('然');
title('[Ca^{2+}]');
cb = colorbar;
ylabel(cb,'on rate');
colormap(rMap);
caxis(onRate([1 end]));
set(gca,'ColorScale','log');
%legend(cellfun(@(c) sprintf('Kon: 1x10^{%d} /M/s',c),num2cell(log10(onRate)),'uni',0),'location','northeast');
set(gca,'fontsize',16);

subplot(1,3,2);
plot(1000*tvec, 1e6*squeeze(values(:,2,plotIdx)), 'linewidth',1.5);
xlim([0 200]);
yLimBuffer = ylim;
xlabel('Time (ms)');
ylabel('然');
title('[CaOGB]');
set(gca,'fontsize',16);

axes('Position',[0.5 0.6 0.1 0.25]);
box on;
plot(1000*tvec, 1e6*squeeze(values(:,2,plotIdx)), 'linewidth',1.5);
xlim([0 10]);
%ylim([0 yLIM(2)]);
xlabel('Time (ms)');
ylabel('然');
set(gca,'fontsize',12);


subplot(1,3,3);
plot(onRate, 1e6*pkCalcium, 'linewidth',1,'marker','o','markerfacecolor','k','color','k');
set(gca,'xscale','log');
xlim(onRate([1 end]));
xlabel('On Rate (1/M/s)');
ylabel('[Ca^{2+}] (然)');
title('K_{on} vs. \Delta[Ca^{2+}]');
set(gca,'fontsize',16);

% This is the fluorophore + instant buffer plots
% ODE Parameters
tspan = [0,0.2];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.02);
dt = 0.00001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 1e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 70e-9; % 70nM
amplitude = 7e-12; % 7pA

% Properties of OGB
numOnRate = 75;
onRateRange = [4 12];
onRate = logspace(onRateRange(1),onRateRange(2),numOnRate); % 1/M/s
ogbConcentration = 50e-6;% 50然 example
ogbKD = 205e-9; % 205nM

fluorProperties.totalConcentration = ogbConcentration;
fluorProperties.kd = ogbKD;

fluorKappa = ogbConcentration / (restConcentration + ogbKD);
initialState = restConcentration * [1 fluorKappa];

sols = cell(1,numOnRate);
values = zeros(NT,2,numOnRate); 
for i = 1:numOnRate
    fluorProperties.onRate = onRate(i);
    
    prmODE = @(t,y) instantPlusFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,20,fluorProperties);
    sols{i} = ode23s(prmODE, tspan, initialState, odeOptions);
    
    values(:,:,i) = deval(sols{i},tvec)';
end

pkCalcium = max(squeeze(values(:,1,:)),[],1); % Max of calcium transient

plotIdx = 1:3:numOnRate;
rMap = [linspace(0.2,1,length(plotIdx))',zeros(length(plotIdx),2)];
figure(2); clf;
set(gcf,'units','normalized','outerposition',[0.2 0.1 0.7 0.4]);
set(gcf,'DefaultAxesColorOrder',rMap);

subplot(1,3,1);
plot(1000*tvec, 1e6*squeeze(values(:,1,plotIdx)), 'linewidth',1.5);
xlim([0 25]);
% ylim(yLimCalcium);
xlabel('Time (ms)');
ylabel('然');
title('[Ca^{2+}]');
cb = colorbar;
ylabel(cb,'on rate');
colormap(rMap);
caxis(onRate([1 end]));
set(gca,'ColorScale','log');
%legend(cellfun(@(c) sprintf('Kon: 1x10^{%d} /M/s',c),num2cell(log10(onRate)),'uni',0),'location','northeast');
set(gca,'fontsize',16);

subplot(1,3,2);
plot(1000*tvec, 1e6*squeeze(values(:,2,plotIdx)), 'linewidth',1.5);
xlim([0 200]);
ylim(yLimBuffer);
xlabel('Time (ms)');
ylabel('然');
title('[CaOGB]');
set(gca,'fontsize',16);

axes('Position',[0.5 0.6 0.1 0.25]);
box on;
plot(1000*tvec, 1e6*squeeze(values(:,2,plotIdx)), 'linewidth',1.5);
xlim([0 10]);
ylim(yLimBuffer);
xlabel('Time (ms)');
ylabel('然');
set(gca,'fontsize',12);

subplot(1,3,3);
plot(onRate, 1e6*pkCalcium, 'linewidth',1,'marker','o','markerfacecolor','k','color','k');
set(gca,'xscale','log');
xlim(onRate([1 end]));
xlabel('On Rate (1/M/s)');
ylabel('[Ca^{2+}] (然)');
title('K_{on} vs. \Delta[Ca^{2+}]');
set(gca,'fontsize',16);





%% ODE for 1 spine with fluorescent buffer only
function dydt = justFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,fluorProperties)
    freeCalcium = y(1); 
    boundBuffer = y(2); 
    freeBuffer = fluorProperties.totalConcentration - boundBuffer;
    
    % Calcium current stuff
    restCurrent = extrusionRate * restConcentration;
    curr = restCurrent + ica(t,amplitude)/(2*96485*spineVolume);
    extrusion = extrusionRate * freeCalcium; 

    % Fluorescent Buffer Reaction
    fluorAssociation = freeCalcium * freeBuffer * fluorProperties.onRate; 
    fluorDissociation = boundBuffer * (fluorProperties.onRate * fluorProperties.kd);
    
    %Output
    dydt(1,1) = curr - extrusion - fluorAssociation + fluorDissociation;
    dydt(2,1) = fluorAssociation - fluorDissociation;
end


%% ODE for 1 spine with fluorescent buffer and instant buffer
function dydt = instantPlusFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,kappa,fluorProperties)
    freeCalcium = y(1); 
    boundBuffer = y(2); 
    freeBuffer = fluorProperties.totalConcentration - boundBuffer;
    
    % Calcium current stuff
    restCurrent = extrusionRate * restConcentration;
    curr = restCurrent + ica(t,amplitude)/(2*96485*spineVolume);
    extrusion = extrusionRate * freeCalcium; 

    % Fluorescent Buffer Reaction
    fluorAssociation = freeCalcium * freeBuffer * fluorProperties.onRate; 
    fluorDissociation = boundBuffer * (fluorProperties.onRate * fluorProperties.kd);
    
    %Output
    dydt(1,1) = (curr - extrusion - fluorAssociation + fluorDissociation) / (1 + kappa);
    dydt(2,1) = fluorAssociation - fluorDissociation;
end









































































