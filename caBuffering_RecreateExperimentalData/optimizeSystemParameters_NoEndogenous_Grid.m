
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
odePrms.tspan = [0,0.4];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.02);
dt = 0.00001;
tvec = odePrms.tspan(1):dt:odePrms.tspan(2);
NT = length(tvec);

% fluorescent buffers
fluorNames = {'fluo5f','ogb1-20然','ogb1-50然','ogb1-100然','fluo4-20然'};
useBufferIdx = logical([0 1 1 1 1]);
% fluorescent parameters
fluorKon = 5e10;
fluorKd = [800e-9, 205e-9, 205e-9, 205e-9, 300e-9];
fluorConcentrations = [300e-6, 20e-6, 50e-6, 100e-6, 20e-6];
fluorPeak = [77e-9, 400e-9, 162.7e-9, 68.9e-9, 474.5e-9];
fluorTau = [100, 30.7, 67.5, 141.1, 20.77]; 

% fluorescent parameter structure
fluorPrms.names = fluorNames;
fluorPrms.useBufferIdx = useBufferIdx;
fluorPrms.kd = fluorKd;
fluorPrms.kon = fluorKon;
fluorPrms.bt = fluorConcentrations;
fluorPrms.peak = fluorPeak;
fluorPrms.tau = fluorTau;

% fit parameters for decay - 
fitPrms.startFitOffset = 0.005; % 5ms
fitPrms.pkEndProportion = 0.1;

% optimizationParameters - Range
NRest = 8;
NBeta = 4;
NAmp = 4;
rest = linspace(40e-9,125e-9,NRest);
beta = linspace(200,3000,NBeta);
iAmp = logspace(-12,-10,NAmp);

costs = zeros(NRest,NBeta,NAmp);
estimates = cell(NRest,NBeta,NAmp);
parfor nbeta = 1:NBeta
    for nrest = 1:NRest
        fprintf('Parfor Loop (%d) on nrest: %d/%d\n',nbeta,nrest,NRest);
        for namp = 1:NAmp
            currentOptPrms = [rest(nrest) beta(nbeta) iAmp(namp) kon(nkon) kd(nkd) bt(nbt)]; %#ok - have to
            [currentCost,estimates{nrest,nbeta,namp}] = ...
                optimizeSystem(currentOptPrms,systemPrms,fluorPrms,odePrms,fitPrms,odeOptions);
            costs(nrest,nbeta,namp) = sum(currentCost.^2);
        end
    end
end
% save('optimizeParametersGridSearch_190815');


%% Report Best Parameters:
currentOptPrms = struct();

[~,costIdx] = sort(costs(:));
[restIdx,betaIdx,ampIdx] = ind2sub(size(costs),costIdx);
bestPrms = [rest(restIdx(1)) beta(betaIdx(1)) iAmp(ampIdx(1))];

NO = sum(fluorPrms.useBufferIdx);
objIdx = find(fluorPrms.useBufferIdx);
bestSols = zeros(NT,2,NO);
for no = 1:NO
    cObjIdx = objIdx(no);
    currentOptPrms.rest= bestPrms(1);
    currentOptPrms.beta = bestPrms(2);
    currentOptPrms.iCalciumAmp = bestPrms(3);

    currentFluorPrms.kon = fluorPrms.kon; % For now this is a single value
    currentFluorPrms.kd = fluorPrms.kd(cObjIdx);
    currentFluorPrms.bt = fluorPrms.bt(cObjIdx);

    cKappaFB = currentFluorPrms.bt / (bestPrms(1) + currentFluorPrms.kd);
    cu0 = [bestPrms(1) bestPrms(1)*cKappaFB];

    odeFunction = @(t,y) justFluorophores(t,y,systemPrms,currentFluorPrms,currentOptPrms);
    cSol = ode23s(odeFunction, odePrms.tspan, cu0, odeOptions);
    bestSols(:,:,no) = deval(cSol,tvec)';
end

% Make Report Legend
reportLegend = cell(NO,1);
for no = 1:NO
    deltaFluor = bestSols(:,2,no) - bestSols(1,2,no);
    [fluorMax,fluorPkIdx] = max(deltaFluor);
    fitStartIdx = find(tvec>= (tvec(fluorPkIdx)+fitPrms.startFitOffset),1,'first'); % start a little after peak
    fluorPkEndValue = fluorMax*fitPrms.pkEndProportion; % Value at which deltaFluor returned to some % of peak
    fitEndIdx = find(deltaFluor > fluorPkEndValue,1,'last'); % Idx at some % of peak

    dx = tvec(fitStartIdx+1:fitEndIdx) - tvec(fitStartIdx);
    fStartPk = deltaFluor(fitStartIdx);
    estTau = 1000*mean(dx' ./ log(fStartPk./deltaFluor(fitStartIdx+1:fitEndIdx)));
    
    reportLegend{no} = sprintf('tau:%.2fms, est:%.2fms',fluorTau(objIdx(no)),estTau);
end


f = figure(170); clf;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.7 0.4]);
set(gcf,'DefaultAxesColorOrder',cat(2, linspace(0,1,NO)',zeros(NO,1),linspace(1,0,NO)'));
xLimMax = [50 400 400];

names = {'Calcium','EndogenousBuffer','FluorescentBuffer'};
legendNames = fluorNames(fluorPrms.useBufferIdx);
for i = 1:2
    subplot(1,2,i); hold on;
    plot(1000*tvec, 1e6*squeeze(bestSols(:,i,:)), 'linewidth', 1.5);
    xlim([0 xLimMax(i)]);
    xlabel('Time (ms)'); 
    ylabel('[X] (然)');
    title(names{i});
    if i==1
        legend(legendNames,'location','northeast'); 
    end
    set(gca,'fontsize',16);
    
    if i==1
        bestReport = {sprintf('Rest: %dnM',bestPrms(1)*1e9);
            sprintf('Beta: %d/s',bestPrms(2));
            sprintf('iAmp: %dpA',bestPrms(3)*1e12)};
        XLIM = xlim;
        YLIM = ylim;
        text(mean(XLIM),mean(YLIM)+diff(YLIM)/5,bestReport);
    end
    
    if i==2
        %kappaFB = fluorPrms.bt(objIdx) ./ (bestPrms(1) + fluorPrms.kd(objIdx));
        %for no = 1:NO
        %    line([0 xLimMax(i)],bestSols(1,3,no)+fluorPrms.amp
        %end
        legend(reportLegend,'location','best');
    end
end



%% ODEs used above

function [cost,estimate] = optimizeSystem(optPrms,systemPrms,fluorPrms,odePrms,fitPrms,odeOptions)
    NO = sum(fluorPrms.useBufferIdx); % How many objectives are we running
    objIdx = find(fluorPrms.useBufferIdx);
    
    referenceValues = cat(1, fluorPrms.peak(fluorPrms.useBufferIdx)', fluorPrms.tau(fluorPrms.useBufferIdx)');
    estimateValues = zeros(NO,2);
    
    % For each objective (different fluorophores)
    for no = 1:NO
        cObjIdx = objIdx(no);
        
        % Make the current optimization (grid) parameter structure
        currentOptPrms.rest= optPrms(1);
        currentOptPrms.beta = optPrms(2);
        currentOptPrms.iCalciumAmp = optPrms(3);
        
        % Load in the current fluorophore
        currentFluorPrms.kon = fluorPrms.kon; % For now this is a single value
        currentFluorPrms.kd = fluorPrms.kd(cObjIdx);
        currentFluorPrms.bt = fluorPrms.bt(cObjIdx);
        
        % Compute kappas for initial state at the current [Ca]rest value
        cKappaFB = currentFluorPrms.bt / (optPrms(1) + currentFluorPrms.kd);
        cu0 = [optPrms(1) optPrms(1)*cKappaFB];
        
        % Parameterize ODE function and compute it
        odeFunction = @(t,y) justFluorophores(t,y,systemPrms,currentFluorPrms,currentOptPrms);
        [cTime, cVals] = ode23s(odeFunction, odePrms.tspan, cu0, odeOptions);
        
        % Estimate Peak (from [CaB]/[Bt] as would be done in an experiment)
        % fractionBound * Kd / (1 - fractionBound) = Calcium --- (at eq.)
        fractionRest = cVals(1,2)/currentFluorPrms.bt;
        fractionPeak = max(cVals(:,2))/currentFluorPrms.bt;
        estimateAtRest = fractionRest*currentFluorPrms.kd / (1 - fractionRest);
        estimateAtPeak = fractionPeak*currentFluorPrms.kd / (1 - fractionPeak);
        estimateValues(no,1) = estimateAtPeak - estimateAtRest; % delta calcium at peak
        
        % Estimate Time Constant from decay
        deltaFluor = cVals(:,2) - cVals(1,2); % delta [CaB] over time
        [fluorMax,fluorPkIdx] = max(deltaFluor); % max, max index
        fitStartIdx = find(cTime>= (cTime(fluorPkIdx)+fitPrms.startFitOffset),1,'first'); % start a little after peak
        fluorPkEndValue = fluorMax*fitPrms.pkEndProportion; % Value at which deltaFluor returned to some % of peak
        fitEndIdx = find(deltaFluor > fluorPkEndValue,1,'last'); % Idx at some % of peak
        
        dx = cTime(fitStartIdx+1:fitEndIdx) - cTime(fitStartIdx);
        fStartPk = deltaFluor(fitStartIdx);
        estimateValues(no,2) = 1000*mean(dx ./ log(fStartPk./deltaFluor(fitStartIdx+1:fitEndIdx)));
    end
    cost = (estimateValues(:) - referenceValues) ./ referenceValues;
    estimate = estimateValues;
end

% ODE for endogenous + fluorescent buffer and optimization of system prms
function dydt = justFluorophores(t,y,systemPrms,fluorPrms,optPrms)
    % Calcium current stuff
    restCurrent = optPrms.beta * optPrms.rest;
    curr = restCurrent + ica(t,optPrms.iCalciumAmp)/(2*96485*systemPrms.v);
    extrusion = optPrms.beta * y(1);

    % Fluorescent Buffer Reaction
    fluorAssociation = y(1)*(fluorPrms.bt - y(2))*fluorPrms.kon;
    fluorDissociation = y(2)*fluorPrms.kon*fluorPrms.kd;
    fluorExchange = fluorAssociation - fluorDissociation;

    %Output
    dydt(1,1) = curr - extrusion - fluorExchange;
    dydt(2,1) = fluorExchange;
end


% ODE for instant, never-saturated buffer
function dy = instantBufferDynamics(t,y,prm,iCalciumAmp)
    % Current/Extrusion Term
    restCurrent = prm.beta/prm.ks * prm.rest;
    curr = restCurrent + ica(t,iCalciumAmp)/(2*96485*prm.v*prm.ks);
    extrusion = prm.beta/prm.ks * y(1);
    % Derivative
    dy = curr - extrusion;
end












































































