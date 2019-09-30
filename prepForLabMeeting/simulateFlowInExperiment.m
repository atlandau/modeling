
hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/prepForLabMeeting';
lmPath = fullfile(cdSL,'presentations/LabMeeting/190906');
dacPath = '/Users/LandauLand/Documents/Research/SabatiniLab/presentations/DACs/DAC2';

twoPanelSize = [0.29 0.27 0.21 0.38];
threePanelSize = [0.29 0.27 0.21*2/3 0.38*2/3];


%% -- Kon too low for instant to be true

tspan = [0,0.5];
odeOptions = odeset('RelTol',1e-12,'AbsTol',1e-12);%,'MaxStep',0.02);
dt = 0.00001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 0.5e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 50e-9; % 50nM
amplitude = 2e-12; % 7pA
sysProperties = [spineVolume, extrusionRate, restConcentration, amplitude];

kappaFunc = @(prop,ca) prop.totalConcentration / (ca + prop.kd);

% Fluorescent Buffer Properties -- fluo-5f
fluo5f.onRate = 5e8;
fluo5f.totalConcentration = 300e-6;
fluo5f.kd = 2.3e-6;
fluo4ff.onRate = 5e8;
fluo4ff.totalConcentration = 600e-6;
fluo4ff.kd = 9e-6;

fluorProperties = fluo4ff;
fluorKappa = kappaFunc(fluorProperties,restConcentration);

% Endogenous Buffer Properties
bufferProperties.onRate = 5e8;
bufferProperties.totalConcentration = 20e-6;
bufferProperties.kd = bufferProperties.totalConcentration/20;


% Vary Fluorophore
varyFluorophore = 40e-6:80e-6:fluorProperties.totalConcentration;
NF = length(varyFluorophore);


% Fit Parameters
ftStart = 10e-3;
ftStartIdx = find(tvec>=ftStart,1);
ftTime = tvec(ftStartIdx:end) - tvec(ftStartIdx(1));
ft = fittype('a*exp(-x/tau)');
fo = fitoptions(ft);
fo.StartPoint = [1 0.01];
fo.Lower = [0 0];
fo.Upper = [2 1];
fo.TolFun = 1e-12;
fo.TolX = 1e-12;

msg = '';
values = zeros(NT,3,NF);
peak = nan(2,NF); % true calcium and estimate
tau = nan(1,NF); % decay from fluorophore
eKappa = nan(1,NF); % estimate of kappa (true estimate from rest)
for nf = 1:NF
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf('Fluorophore Concentration %d/%d ... \n',nf,NF);
    fprintf(1,msg);
    fluorProperties.totalConcentration = varyFluorophore(nf);
    fluorKappa = kappaFunc(fluorProperties,restConcentration);
    iState = restConcentration * [1 kappaFunc(bufferProperties,restConcentration) fluorKappa];
            
    prmODE = @(t,y) endogenousAndFluorescent(t,y,sysProperties,bufferProperties,fluorProperties);
    csol = ode23s(prmODE,tspan,iState,odeOptions);
    values(:,:,nf) = deval(csol,tvec)';
    
    if fluorProperties.totalConcentration == 0, continue, end
    
    % Estimate Peak / Tau
    peak(1,nf) = max(values(:,1,nf)) - restConcentration;
    fOccupiedRest = values(1,3,nf) / fluorProperties.totalConcentration;
    fOccupiedPeak = max(values(:,3,nf))/fluorProperties.totalConcentration;
    estimateRest = fOccupiedRest*fluorProperties.kd / (1 - fOccupiedRest);
    estimatePeak = fOccupiedPeak*fluorProperties.kd / (1 - fOccupiedPeak);
    peak(2,nf) = estimatePeak - estimateRest;
    
    cfit = fit(ftTime(:),norman(values(ftStartIdx:end,3,nf)),ft,fo);
    tau(nf) = cfit.tau;
    
    eKappa(nf) = fluorKappa;
end


apCoeff = [ones(NF,1) eKappa(:)] \ (1./peak(2,:))';
tauCoeff = [ones(NF,1) eKappa(:)] \ tau(:);
fprintf(1,'1/AP estimate: %.1f\n',apCoeff(1)/apCoeff(2));
fprintf(1,'tau estimate: %.1f\n',tauCoeff(1)/tauCoeff(2));


figure(1); clf;

subplot(1,2,1); hold on;
plot(eKappa, 1./peak(2,:),'color','k','marker','*','linewidth',1);
line(eKappa([1 end]), apCoeff(2)*eKappa([1 end])+apCoeff(1), 'color','r','linewidth',1);
xlabel('Kappa (from rest)');
ylabel('1/AP');
set(gca,'fontsize',18);

subplot(1,2,2); hold on;
plot(eKappa, tau,'color','k','marker','*','linewidth',1);
line(eKappa([1 end]), tauCoeff(2)*eKappa([1 end])+tauCoeff(1), 'color','r','linewidth',1);
xlabel('Kappa (from rest)');
ylabel('Tau');
set(gca,'fontsize',18);






%% Vary Kon, Bt/Kd, and Kd to lower and lower values

tspan = [0,0.3];
odeOptions = odeset('RelTol',1e-12,'AbsTol',1e-12);%,'MaxStep',0.02);
dt = 0.00005;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 0.5e-15; % 1fL
extrusionRate = 1500; % 1/s 
restConcentration = 50e-9; % 50nM
amplitude = 2e-12; % 7pA
sysProperties = [spineVolume, extrusionRate, restConcentration, amplitude];

kappaFunc = @(prop,ca) prop.totalConcentration / (ca + prop.kd);

% Fluorescent Buffer Properties -- fluo-5f
fluo5f.onRate = 5e8;
fluo5f.totalConcentration = 300e-6;
fluo5f.kd = 2.3e-6;
fluo4ff.onRate = 5e8;
fluo4ff.totalConcentration = 600e-6;
fluo4ff.kd = 9e-6;

fluorProperties = fluo5f;
fluorKappa = kappaFunc(fluorProperties,restConcentration);

% Vary Fluorophore
varyFluorophore = 20e-6:20e-6:fluorProperties.totalConcentration;
NF = length(varyFluorophore);

% Endogenous Buffer Properties
numValues = 15;
onRateRange = [4 9];
btRange = [-7 -2]; % note that this varies Bt and Kd together
kdRange = [-8 -3];
onRate = logspace(onRateRange(1),onRateRange(2),numValues);
bt = logspace(btRange(1),btRange(2),numValues);
kd = logspace(kdRange(1),kdRange(2),numValues);

% Base Endogenous Buffer
bufferBase.onRate = 5e8;
bufferBase.totalConcentration = 4e-3;
bufferBase.kd = bufferBase.totalConcentration/20;
bufferProperties = bufferBase;


% Fit Parameters
ftStart = 10e-3;
ftStartIdx = find(tvec>=ftStart,1);
ftTime = tvec(ftStartIdx:end) - tvec(ftStartIdx(1));
ft = fittype('a*exp(-x/tau)');
fo = fitoptions(ft);
fo.StartPoint = [1 0.01];
fo.Lower = [0 0];
fo.Upper = [2 1];
fo.TolFun = 1e-12;
fo.TolX = 1e-12;

performExperiments = [1 2];
NE = length(performExperiments);

msg = '';
valuesBase = zeros(NT,3,NF);
peakBase = nan(2,NF); % true calcium and estimate
tauBase = nan(1,NF); % decay from fluorophore
eKappaBase = nan(1,NF); % estimate of kappa (true estimate from rest)
for nf = 1:NF
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf('Running experiment for optimal buffer: concentration %d/%d ... \n',nf,NF);
    fprintf(1,msg);
    
    fluorProperties.totalConcentration = varyFluorophore(nf);
    fluorKappa = kappaFunc(fluorProperties,restConcentration);
    iState = restConcentration * [1 kappaFunc(bufferProperties,restConcentration) fluorKappa];

    prmODE = @(t,y) endogenousAndFluorescent(t,y,sysProperties,bufferProperties,fluorProperties);
    csol = ode23s(prmODE,tspan,iState,odeOptions);
    valuesBase(:,:,nf) = deval(csol,tvec)';

    if fluorProperties.totalConcentration == 0, continue, end

    % Estimate Peak / Tau
    peakBase(1,nf) = max(valuesBase(:,1,nf)) - restConcentration;
    fOccupiedRest = valuesBase(1,3,nf) / fluorProperties.totalConcentration;
    fOccupiedPeak = max(valuesBase(:,3,nf))/fluorProperties.totalConcentration;
    estimateRest = fOccupiedRest*fluorProperties.kd / (1 - fOccupiedRest);
    estimatePeak = fOccupiedPeak*fluorProperties.kd / (1 - fOccupiedPeak);
    peakBase(2,nf) = estimatePeak - estimateRest;

    cfit = fit(ftTime(:),norman(valuesBase(ftStartIdx:end,3,nf)),ft,fo);
    tauBase(1,nf) = cfit.tau;

    eKappaBase(1,nf) = fluorKappa;
end


msg = '';
values = zeros(NT,3,NF,NE,numValues);
peak = nan(2,NF,NE,numValues); % true calcium and estimate
tau = nan(1,NF,NE,numValues); % decay from fluorophore
eKappa = nan(1,NF,NE,numValues); % estimate of kappa (true estimate from rest)
apCoeff = zeros(2,NE,numValues);
tauCoeff = zeros(2,NE,numValues);
endogenousEstimate = zeros(2,NE,numValues);
for experiment = 1:NE
    bufferProperties = bufferBase;
    for nv = 1:numValues
        fprintf(1,repmat('\b',1,length(msg)));
        msg = sprintf('Experiment %d/%d, Value %d/%d ... \n',experiment,NE,nv,numValues);
        fprintf(1,msg);
        switch performExperiments(experiment)
            case 1
                bufferProperties.onRate = onRate(nv);
            case 2 
                bufferProperties.totalConcentration = bt(nv);
                bufferProperties.kd = bt(nv)/20;
            case 3 
                bufferProperties.kd = kd(nv)/20;
        end
        for nf = 1:NF
            fluorProperties.totalConcentration = varyFluorophore(nf);
            fluorKappa = kappaFunc(fluorProperties,restConcentration);
            iState = restConcentration * [1 kappaFunc(bufferProperties,restConcentration) fluorKappa];

            prmODE = @(t,y) endogenousAndFluorescent(t,y,sysProperties,bufferProperties,fluorProperties);
            csol = ode23s(prmODE,tspan,iState,odeOptions);
            values(:,:,nf,experiment,nv) = deval(csol,tvec)';

            if fluorProperties.totalConcentration == 0, continue, end

            % Estimate Peak / Tau
            peak(1,nf,experiment,nv) = max(values(:,1,nf,experiment,nv)) - restConcentration;
            fOccupiedRest = values(1,3,nf,experiment,nv) / fluorProperties.totalConcentration;
            fOccupiedPeak = max(values(:,3,nf,experiment,nv))/fluorProperties.totalConcentration;
            estimateRest = fOccupiedRest*fluorProperties.kd / (1 - fOccupiedRest);
            estimatePeak = fOccupiedPeak*fluorProperties.kd / (1 - fOccupiedPeak);
            peak(2,nf,experiment,nv) = estimatePeak - estimateRest;

            cfit = fit(ftTime(:),norman(values(ftStartIdx:end,3,nf,experiment,nv)),ft,fo);
            tau(1,nf,experiment,nv) = cfit.tau;

            eKappa(1,nf,experiment,nv) = fluorKappa;
        end
        
        % Analyze
        apCoeff(:,experiment,nv) = [ones(NF,1) eKappa(1,:,experiment,nv)'] \ (1./peak(2,:,experiment,nv)');
        tauCoeff(:,experiment,nv) = [ones(NF,1) eKappa(1,:,experiment,nv)'] \ (tau(1,:,experiment,nv)');
        endogenousEstimate(1,experiment,nv) = apCoeff(1,experiment,nv)/apCoeff(2,experiment,nv);
        endogenousEstimate(2,experiment,nv) = tauCoeff(1,experiment,nv)/tauCoeff(2,experiment,nv);
        
        xLIM = [-25 140];
        
        figure(1); clf;
        set(gcf,'units','normalized','outerposition',[0.10 0.3 0.51 0.68]);
        
        for i = 1:3, subplot(3,3,i); plotvc(tvec, norman(squeeze(values(:,i,:,experiment,nv))),gcf); end
        subplot(3,3,1);
        title(sprintf('Experiment %d, Value %d', experiment, nv));
        
        for i = 1:3, subplot(3,3,i+3); plotvc(tvec, squeeze(values(:,i,:,experiment,nv)),gcf); end
        
        subplot(3,3,7); 
        plot(eKappa(1,:,experiment,nv),1./peak(2,:,experiment,nv),'color','k','marker','*');
        line(xLIM,xLIM*apCoeff(2,experiment,nv)+apCoeff(1,experiment,nv),'color','r','linewidth',0.5);
        line(xLIM,[0 0],'color','k','linewidth',0.5);
        xlim(xLIM);
        title('1/Peak');
        
        subplot(3,3,8); 
        plot(eKappa(1,:,experiment,nv),tau(1,:,experiment,nv),'color','k','marker','*');
        line(xLIM,xLIM*tauCoeff(2,experiment,nv)+tauCoeff(1,experiment,nv),'color','r','linewidth',0.5);
        line(xLIM,[0 0],'color','k','linewidth',0.5);
        xlim(xLIM);
        title('1/Tau');
        
        subplot(3,3,9);
        bar([1 2 3],[20 endogenousEstimate(:,experiment,nv)']);
        set(gca,'xticklabel',{'true','1/peak','tau'});
        
        print(gcf,'-painters',fullfile(cd,'resultsToday',sprintf('exp%d-value%d',experiment,nv)),'-djpeg');
    end
end

%%
cols = 'kbr';
figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.06 0.48 0.38 0.44]);
subplot(1,2,1); hold on;
for experiment = 1:NE
    plot(1:numValues,squeeze(endogenousEstimate(1,experiment,:)),'color',cols(experiment),'linewidth',1.5,'marker','*');
end

subplot(1,2,2); hold on;
for experiment = 1:NE
    plot(1:numValues,squeeze(endogenousEstimate(2,experiment,:)),'color',cols(experiment),'linewidth',1.5,'marker','*');
end

data2plot = values;
data2plot = norman(data2plot);
figure(2); clf;
set(gcf,'units','normalized','outerposition',[0.44 0.48 0.52 0.44]);
for i = 1:3
    subplot(1,3,i); hold on;
    plotvc(tvec, squeeze(data2plot(:,i,:,2,1)),gcf);
    xlim([0 0.1]);
end

%% ODE for 1 spine with endogenous and fluorescent buffer only
function dydt = endogenousAndFluorescent(t,y,sysProperties,bufferProperties,fluorProperties)
    freeCalcium = y(1);
    boundBuffer = y(2);
    boundFluor = y(3);
    freeBuffer = bufferProperties.totalConcentration - boundBuffer;
    freeFluor = fluorProperties.totalConcentration - boundFluor;
    
    % Calcium current stuff
    restCurrent = sysProperties(2) * sysProperties(3);
    curr = restCurrent + ica(t,sysProperties(4))/(2*96485*sysProperties(1));
    extrusion = sysProperties(2) * freeCalcium; 

    % Endogenous Buffer Reaction 
    bufferAssociation = freeCalcium * freeBuffer * bufferProperties.onRate; 
    bufferDissociation = boundBuffer * (bufferProperties.onRate * bufferProperties.kd);
    bufferExchange = bufferAssociation - bufferDissociation;
    
    % Fluorescent Buffer Reaction
    fluorAssociation = freeCalcium * freeFluor * fluorProperties.onRate;
    fluorDissociation = boundFluor * (fluorProperties.onRate * fluorProperties.kd);
    fluorExchange = fluorAssociation - fluorDissociation;
    
    %Output
    dydt(1,1) = curr - extrusion - bufferExchange - fluorExchange;
    dydt(2,1) = bufferExchange;
    dydt(3,1) = fluorExchange;
end


































































