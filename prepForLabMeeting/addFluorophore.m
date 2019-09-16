
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


%% -- vary Bt and Kd to maintain kappa=20 and do a bunch of Kons

% ODE Parameters
tspan = [0,0.02];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',0.02);
dt = 0.0001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 1e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 50e-9; % 50nM
amplitude = 2e-12; % 7pA

% Properties of Endogenous Buffer
numOnRate = 50;
onRateRange = [4 12];
onRate = logspace(onRateRange(1),onRateRange(2),numOnRate); % 1/M/s
% onRate = linspace(10^onRateRange(1),10^onRateRange(2),numOnRate); % 1/M/s

numBT = 10;
kappa = 20;
btRateRange = [-4 -1];
bt = logspace(btRateRange(1),btRateRange(2),numBT);
kd = bt/kappa;

% Get Instant Solution
systemPrms.beta=extrusionRate;
systemPrms.rest=restConcentration;
systemPrms.ks=kappa;
systemPrms.v=spineVolume;
systemPrms.amp = amplitude;

instantSol = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms),tspan,restConcentration,odeOptions);
trueSol = deval(instantSol,tvec)';

msg = '';
fprintf(1,'working... \n');
sols = cell(numOnRate,numBT);
values = zeros(NT,2,numOnRate,numBT);
sse = zeros(numOnRate,numBT);
for non = 1:numOnRate
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf('OnRate %d/%d...\n',non,numOnRate);
    fprintf(1,msg);
    for nbt = 1:numBT
        fluorProperties.totalConcentration = bt(nbt);
        fluorProperties.kd = kd(nbt);
        fluorProperties.onRate = onRate(non);
        initialState = [restConcentration restConcentration*bt(nbt)/(restConcentration+kd(nbt))];
        prmODE = @(t,y) justFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,fluorProperties);
        sols{non,nbt} = ode23s(prmODE,tspan,initialState,odeOptions);
        values(:,:,non,nbt) = deval(sols{non,nbt},tvec)';
        sse(non,nbt) = sum(((values(:,1,non,nbt) - trueSol)./trueSol).^2);
    end
end


%% ----- plot kon traces and sse -----
figure(189); clf;
set(gcf,'units','normalized','outerposition',twoPanelSize);
hold on;
skipNumber = 0;
cmap = [linspace(0.1,1,numOnRate-skipNumber)'.^(1/2),zeros(numOnRate-skipNumber,2)];
for non = 1:1:numOnRate-skipNumber
    plot(tvec,values(:,1,non,end),'color',cmap(non,:),'linewidth',1.5);
end
plot(tvec,trueSol,'color','k','linewidth',3);
set(gca,'xtick',0:0.005:0.02);
set(gca,'xticklabel',1000*get(gca,'xtick'));
set(gca,'ytick',(0:1:3)*1e-6);
set(gca,'yticklabel',1e6*get(gca,'ytick'));
xlabel('Time (ms)');
ylabel('[Ca^{2+}]');
set(gca,'fontsize',18);
cb = colorbar('Location','Manual','Position',[0.7 0.2 0.08 0.68]);
colormap(gca,cmap);
cb.Ticks = [0 1];
cb.TickLabels = fliplr({sprintf('1x10^{%d}',onRateRange(2)),sprintf('1x10^{%d}',onRateRange(1))});

[~,kon5x8] = min(abs(onRate-5e8));
inset = axes('Position',[0.3 0.45 0.35 0.35]);
hold(inset,'on');
plot(tvec,trueSol,'color','k','linewidth',1);
plot(tvec,values(:,1,kon5x8,end),'color',cmap(kon5x8,:),'linewidth',1);
set(inset,'Visible','off');
% print(gcf,'-painters',fullfile(dacPath,'OnRateDemonstration'),'-djpeg');


figure(190); clf;
set(gcf,'units','normalized','outerposition',[0.6 0.27 0.21 0.38]);
hold on;
shadedErrorBar(onRate, mean(sse,2), std(sse,[],2)./sqrt(size(sse,2)), {'color','k','linewidth',1.5},0.3);
set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('On-Rate');
ylabel('Sum Squared Error');
set(gca,'fontsize',18);
% print(gcf,'-painters',fullfile(lmPath,'OnRateError'),'-djpeg');




%% Fix Kon, Vary Kd and Bt


% ODE Parameters
tspan = [0,0.02];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',0.02);
dt = 0.0001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 1e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 50e-9; % 50nM
amplitude = 2e-12; % 7pA

% Properties of Endogenous Buffer
onRate = 5e8; 
numBT = 200;
btRange = [-5 0];
bt = logspace(btRange(1),btRange(2),numBT);
numKD = 200;
kdRange = [-7 -1];
kd = logspace(kdRange(1),kdRange(2),numKD);

% Get Instant Solution
systemPrms.beta=extrusionRate;
systemPrms.rest=restConcentration;
systemPrms.ks=kappa;
systemPrms.v=spineVolume;
systemPrms.amp = amplitude;
fluorProperties.onRate = onRate;

instantSol = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms),tspan,restConcentration,odeOptions);
trueSol = deval(instantSol,tvec)';

msg = '';
fprintf(1,'working... \n');
sols = cell(numBT,numKD);
values = zeros(NT,2,numBT,numKD);
sse = zeros(numBT,numKD);
for nbt = 1:numBT
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf('Bt %d/%d...\n',nbt,numBT);
    fprintf(1,msg);
    for nkd = 1:numKD
        fluorProperties.totalConcentration = bt(nbt);
        fluorProperties.kd = kd(nkd);
        initialState = [restConcentration restConcentration*bt(nbt)/(restConcentration+kd(nkd))];
        prmODE = @(t,y) justFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,fluorProperties);
        sols{nbt,nkd} = ode23s(prmODE,tspan,initialState,odeOptions);
        values(:,:,nbt,nkd) = deval(sols{nbt,nkd},tvec)';
        sse(nbt,nkd) = sum(((values(:,1,nbt,nkd) - trueSol)./trueSol).^2);
    end
end

% save('kdVsBt_results190905','sols','values','sse','trueSol','tvec','kd','bt','onRate','systemPrms');

%% Plot beautiful Kd vs BT - SSE plot
figure(202); clf;
set(gcf,'units','normalized','outerposition',[0.60 0.27 0.22 0.38]);
hold on;
h = pcolor(kd,bt,sse);
h.EdgeColor = 'none';
xlim(kd([1 end]));
ylim(bt([1 end]));
set(gca,'xscale','log');
set(gca,'yscale','log');
xTick = get(gca,'xtick');
set(gca,'xtick',xTick(2:2:end));
set(gca,'ytick',get(gca,'ytick'));
set(gca,'xticklabel',readableMolar(get(gca,'xtick')));
set(gca,'yticklabel',readableMolar(get(gca,'ytick')));
xlabel('K_D');
ylabel('[B]_T');
colormap(flipud(hot));
colorbar;
set(gca,'colorscale','log');
set(gca,'fontsize',18);
% print(gcf,'-painters',fullfile(lmPath,'kdVsBt_beautifulFuckingPlot'),'-djpeg');
% plot(kd([1 end]),20*kd([1 end]),'color','b','linewidth',2);



%% Plot SSE vs. Bt for optimal Kd

kappaFit = nan(numBT,1);
bestKD = nan(numBT,1);
for nbt = 1:numBT
    diffBtOptimal = kd - bt(nbt)/kappa;
    if sum(diffBtOptimal<0)==0, continue, end
    if sum(diffBtOptimal>0)==0, continue, end 
    [~,idxKD] = min(abs(diffBtOptimal));
    kappaFit(nbt) = sse(nbt,idxKD);
    bestKD(nbt) = kd(idxKD);
end

figure(203); clf;
set(gcf,'units','normalized','outerposition',[0.60 0.27 0.22 0.38]);

skip=6;
idxNan = find(~isnan(kappaFit));
plot(bt(idxNan(1:skip:end)),kappaFit(idxNan(1:skip:end)),'color','k','linewidth',1.5);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xtick',get(gca,'xtick'));
set(gca,'xticklabel',readableMolar(get(gca,'xtick')));
xlabel('B_T');
ylabel('SSE');
set(gca,'fontsize',18);
% print(gcf,'-painters',fullfile(lmPath,'kdVsBt_bestKappaLineSSEs'),'-djpeg');



%% --------- now do it with scaled kappa -----------

% ODE Parameters
tspan = [0,0.02];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',0.02);
dt = 0.0001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 1e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 50e-9; % 50nM
amplitude = 2e-12; % 7pA

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

% --- turn off fluorophore ---
fluorProperties.totalConcentration = 0;

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
rMap = [linspace(0.2,1,length(plotIdx))',zeros(length(plotIdx),2)].^(1/1.2);

figure(111); clf;
set(gcf,'units','normalized','outerposition',[0.38 0.27 0.22 0.38]);
hold on;
for i = 1:length(plotIdx)
    plot(1000*tvec, 1e6*squeeze(values(:,1,plotIdx(i))), 'linewidth',1.5,'color',rMap(i,:));
end
xlim([0 20]);
yLimCalcium = ylim;
xlabel('Time (ms)');
ylabel('[Ca^{2+}] (然)');
cb = colorbar('Location','Manual','Position',[0.7 0.2 0.08 0.68]);
% ylabel(cb,'kappaKD');
colormap(rMap);
caxis(kappaKD([1 end]));
set(gca,'ytick',0:1:4);
set(gca,'ColorScale','log');
cb.Ticks = cb.Ticks;
cb.TickLabels = readableMolar(cb.Ticks);
set(gca,'fontsize',18);
% print(gcf,'-painters',fullfile(lmPath,'KappaKD_Traces'),'-djpeg');

figure(112); clf;
set(gcf,'units','normalized','outerposition',[0.60 0.27 0.22 0.38]);
hold on;
plot(kappaKD, 1e6*pkCalcium, 'linewidth',1,'marker','o','color','k');
% plot(kappaKD, 1e6*pkEstBuffer, 'linewidth',1,'marker','o','color','r');
set(gca,'xscale','log');
set(gca,'xtick',[1e-8 1e-6 1e-4]);
set(gca,'xticklabel',readableMolar(get(gca,'xtick')));
set(gca,'ytick',0:1:4);
xlim(kappaKD([1 end]));
yLIM = ylim;
xlabel('Kappa Kd');
ylabel('[Ca^{2+}] (然)');
% title('kappaKD vs. \Delta[Ca^{2+}]');
% legend('TruePeak','Estimate','location','northeast');
set(gca,'fontsize',18);
% print(gcf,'-painters',fullfile(lmPath,'KappaKD_PkCalcium'),'-djpeg');




%% --------- scaled kappa plus fluorophore -----------

% ODE Parameters
tspan = [0,0.25];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',0.02);
dt = 0.0001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 1e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 50e-9; % 50nM
amplitude = 2e-12; % 7pA

%numKappaKD = 100;
%kappaRangeKD = [-8 -3];
%kappaKD = logspace(kappaRangeKD(1),kappaRangeKD(2),numKappaKD);
kappaKD = 175e-6;

% On Rate of Fluorophore
onRate = 5e8; 

% Properties of Fluorophores
ogbConcentration = 50e-6;% 50然 example
ogbKD = 205e-9; % 205nM
fluo5fConcentration = 300e-6; 
fluo5fKD = 2.3e-6; %800e-9;
fluo4ffConcentration = 100e-6;
fluo4ffKD = 9e-6;

useFluor = 'fluo5f';
fluorProperties.onRate = onRate;
fluorProperties.totalConcentration = eval([useFluor,'Concentration']);
fluorProperties.kd = eval([useFluor,'KD']);

delayTimeIdx = find(tvec>=0.005,1);

concentrations = (0:25:300)*1e-6;
numConcentrations = length(concentrations);
sols = cell(1,numConcentrations);
values = zeros(NT,2,numConcentrations); 
pkEstBuffer = zeros(numConcentrations,1);
pkEstBuffDelay = zeros(numConcentrations,1);
for i = 1:numConcentrations
    fluorProperties.totalConcentration = concentrations(i);
    fluorKappa = fluorProperties.totalConcentration / (restConcentration + fluorProperties.kd);
    initialState = restConcentration * [1 fluorKappa];
    
    prmODE = @(t,y) instantScaledPlusFluorophores(t,y,spineVolume,extrusionRate,restConcentration,amplitude,kappa,kappaKD,fluorProperties);
    sols{i} = ode23s(prmODE, tspan, initialState, odeOptions);
    
    values(:,:,i) = deval(sols{i},tvec)';
    fOccupiedRest = values(1,2,i)/fluorProperties.totalConcentration;
    fOccupiedPeak = max(values(:,2,i))/fluorProperties.totalConcentration;
    estimateRest = fOccupiedRest*fluorProperties.kd / (1 - fOccupiedRest);
    estimatePeak = fOccupiedPeak*fluorProperties.kd / (1 - fOccupiedPeak);
    pkEstBuffer(i) = estimatePeak - estimateRest;
    
    fOccupiedPkDelay = values(delayTimeIdx,2,i)/fluorProperties.totalConcentration;
    estimatePkDelay = fOccupiedPkDelay*fluorProperties.kd / (1 - fOccupiedPkDelay)';
    pkEstBuffDelay(i) = estimatePkDelay - estimateRest;
end

pkCalcium = max(squeeze(values(:,1,:)),[],1) - restConcentration; % Max of calcium transient
pkDelayCalcium = squeeze(values(delayTimeIdx,1,:))-restConcentration; % Calcium transient value at 5ms (after escape)
caBuffer = squeeze(values(:,2,:));
caBuffer = (caBuffer - repmat(caBuffer(1,:),size(caBuffer,1),1))./repmat(caBuffer(1,:),size(caBuffer,1),1);

plotIdx = 1:numConcentrations;
rMap = [linspace(0.2,1,length(plotIdx))',zeros(length(plotIdx),2)].^(1/1.2);

figure(111); clf;
set(gcf,'units','normalized','outerposition',[0.26 0.27 0.22 0.38]);
hold on;
for i = 1:length(plotIdx)
    plot(1000*tvec, 1e6*squeeze(values(:,1,plotIdx(i))), 'linewidth',1.5,'color',rMap(i,:));
end
xlim([-20 tspan(2)*1000]);
yLimCalcium = ylim;
xlabel('Time (ms)');
ylabel('[Ca^{2+}] (然)');
cb = colorbar('Location','Manual','Position',[0.7 0.2 0.08 0.68]);
% ylabel(cb,'kappaKD');
colormap(rMap);
caxis(concentrations([1 end]));
% set(gca,'ytick',0:1:4);
% set(gca,'ColorScale','log');
cb.Ticks = cb.Ticks;
cb.TickLabels = readableMolar(cb.Ticks);
set(gca,'fontsize',18);

figure(112); clf;
set(gcf,'units','normalized','outerposition',[0.48 0.27 0.22 0.38]);
hold on;
plot(concentrations, 1e6*pkCalcium, 'linewidth',1,'marker','o','color','k');
plot(concentrations, 1e6*pkEstBuffer,'linewidth',1,'marker','o','color','r');
% plot(kappaKD, 1e6*pkEstBuffer, 'linewidth',1,'marker','o','color','r');
% set(gca,'xscale','log');
% set(gca,'xtick',[1e-8 1e-6 1e-4]);
set(gca,'xticklabel',readableMolar(get(gca,'xtick')));
% set(gca,'ytick',0:1:4);
xlim(concentrations([1 end]));
yLIM = ylim;
xlabel('Concentration Fluo-5f');
ylabel('[Ca^{2+}] (然)');
% title('kappaKD vs. \Delta[Ca^{2+}]');
% legend('TruePeak','Estimate','location','northeast');
set(gca,'fontsize',18);


figure(113); clf;
set(gcf,'units','normalized','outerposition',[0.04 0.27 0.22 0.38]);
hold on;
for i = 1:length(plotIdx)
    plot(1000*tvec, 1e6*caBuffer(:,i), 'linewidth',1.5,'color',rMap(i,:));
end
xlim([-20 tspan(2)*1000]);
yLimCalcium = ylim;
xlabel('Time (ms)');
ylabel('[CaFluo5f] (然)');
set(gca,'fontsize',18);

figure(114); clf;
set(gcf,'units','normalized','outerposition',[0.70 0.27 0.22 0.38]);
hold on;
plot(concentrations, 1e6*pkDelayCalcium, 'linewidth',1,'marker','o','color','k');
plot(concentrations, 1e6*pkEstBuffDelay,'linewidth',1,'marker','o','color','r');
% plot(kappaKD, 1e6*pkEstBuffer, 'linewidth',1,'marker','o','color','r');
% set(gca,'xscale','log');
% set(gca,'xtick',[1e-8 1e-6 1e-4]);
set(gca,'xticklabel',readableMolar(get(gca,'xtick')));
% set(gca,'ytick',0:1:4);
xlim(concentrations([1 end]));
yLIM = ylim;
xlabel('Concentration Fluo-5f');
ylabel('[Ca^{2+}] (然)');
% title('kappaKD vs. \Delta[Ca^{2+}]');
% legend('TruePeak','Estimate','location','northeast');
set(gca,'fontsize',18);



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
    dydt(1,1) = (curr - extrusion + fluorDissociation - fluorAssociation)/(1+currentKappa);
    dydt(2,1) = fluorAssociation - fluorDissociation;
end



































































