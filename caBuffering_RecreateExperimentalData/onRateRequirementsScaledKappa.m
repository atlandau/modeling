
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

u0 = systemPrms.rest;
tspan = [0,0.5];
dt = 0.000001;
tvec = tspan(1):dt:tspan(2);

% Instant Solution
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',dt);
[trueTime,trueSol] = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms,iCalciumAmp),tspan,u0,odeOptions);
plot(trueTime,trueSol);


%% Optimize endogenous:[Kon,Kd,Bt], rest, beta, and iCalcium to get following results with fluorophores in the cell

% ODE Parameters
tspan = [0,1];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.02);
dt = 0.00001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 1e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 70e-9; % 70nM
amplitude = 7e-12; % 7pA
kappa = 20;

numKappaKD = 100;
kappaRangeKD = [-8 -3];
kappaKD = logspace(kappaRangeKD(1),kappaRangeKD(2),numKappaKD);

% On Rate of Fluorophore
onRate = 5e8; 

% Properties of Fluorophores
ogbConcentration = 50e-6;% 50然 example
ogbKD = 205e-9; % 205nM
fluo5fConcentration = 100e-6; 
fluo5fKD = 800e-9;
fluo4ffConcentration = 100e-6;
fluo4ffKD = 9e-6;

useFluor = 'fluo5f';
fluorProperties.onRate = onRate;
fluorProperties.totalConcentration = eval([useFluor,'Concentration']);
fluorProperties.kd = eval([useFluor,'KD']);

fluorKappa = fluorProperties.totalConcentration / (restConcentration + fluorProperties.kd);
initialState = restConcentration * [1 fluorKappa];

sols = cell(1,numKappaKD);
values = zeros(NT,2,numKappaKD); 
pkEstBuffer = zeros(numKappaKD,1);
for i = 1:numKappaKD    
    prmODE = @(t,y) instantScaledPlusFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,kappa,kappaKD(i),fluorProperties);
    sols{i} = ode23s(prmODE, tspan, initialState, odeOptions);
    
    values(:,:,i) = deval(sols{i},tvec)';
    fOccupiedRest = values(1,2,i)/fluorProperties.totalConcentration;
    fOccupiedPeak = max(values(:,2,i))/fluorProperties.totalConcentration;
    estimateRest = fOccupiedRest*fluorProperties.kd / (1 - fOccupiedRest);
    estimatePeak = fOccupiedPeak*fluorProperties.kd / (1 - fOccupiedPeak);
    pkEstBuffer(i) = estimatePeak - estimateRest;
end

pkCalcium = max(squeeze(values(:,1,:)),[],1) - restConcentration; % Max of calcium transient

plotIdx = 1:3:numKappaKD;
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
ylabel(cb,'kappaKD');
colormap(rMap);
caxis(kappaKD([1 end]));
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

axes('Position',[0.535 0.65 0.08 0.25]);
box on;
plot(1000*tvec, 1e6*squeeze(values(:,2,plotIdx)), 'linewidth',1.5);
xlim([0 10]);
%ylim([0 yLIM(2)]);
xlabel('Time (ms)');
ylabel('然');
set(gca,'fontsize',12);


subplot(1,3,3); hold on;
plot(kappaKD, 1e6*pkCalcium, 'linewidth',1,'marker','o','color','k');
plot(kappaKD, 1e6*pkEstBuffer, 'linewidth',1,'marker','o','color','r');
set(gca,'xscale','log');
xlim(kappaKD([1 end]));
yLIM = ylim;
% ylim([0 yLIM(2)]);
xlabel('Kappa Kd');
ylabel('[Ca^{2+}] (然)');
title('kappaKD vs. \Delta[Ca^{2+}]');
legend('TruePeak','Estimate','location','northeast');
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
    dydt(1,1) = curr - extrusion - fluorAssociation + fluorDissociation;
    dydt(2,1) = fluorAssociation - fluorDissociation;
    dydt = dydt/(1+kappa); % everything scaled by instant buffer
end

%% ODE for 1 spine with fluorescent buffer and instant buffer that scales kappa
function dydt = instantScaledPlusFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,kappa,kappaKD,fluorProperties)
    freeCalcium = y(1); 
    boundBuffer = y(2); 
    freeBuffer = fluorProperties.totalConcentration - boundBuffer;
    currentKappa = kappa * kappaKD/(freeCalcium + kappaKD);
    
    % Calcium current stuff
    restCurrent = extrusionRate * restConcentration;
    curr = restCurrent + ica(t,amplitude)/(2*96485*spineVolume);
    extrusion = extrusionRate * freeCalcium; 

    % Fluorescent Buffer Reaction
    fluorAssociation = freeCalcium * freeBuffer * fluorProperties.onRate; 
    fluorDissociation = boundBuffer * (fluorProperties.onRate * fluorProperties.kd);
    
    %Output
    dydt(1,1) = (curr - extrusion + fluorDissociation - fluorAssociation)/currentKappa;
    dydt(2,1) = fluorAssociation - fluorDissociation;
end









































































