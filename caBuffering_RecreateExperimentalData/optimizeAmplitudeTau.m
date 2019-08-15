
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
odePrms.tspan = [0,0.5];
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',dt);

% fluorescent buffers
names = {'fluo5f','ogb1','ogb1','ogb1','fluo4'};
useBufferIdx = logical([0 1 1 1 1]);
% fluorescent parameters
fluorKon = 5e8;
fluorKd = [800e-9, 205e-9, 205e-9, 205e-9, 300e-9];
fluorConcentrations = [300e-6, 20e-6, 50e-6, 100e-6, 20e-6];
fluorPeak = [77e-9, 400e-9, 162.7e-9, 68.9e-9, 474.5e-9];
fluorTau = [100, 30.7, 67.5, 141.1, 20.77]; 

% fluorescent parameter structure
fluorPrms.useBufferIdx = useBufferIdx;
fluorPrms.kd = fluorKd;
fluorPrms.kon = fluorKon;
fluorPrms.bt = fluorConcentrations;
fluorPrms.peak = fluorPeak;
fluorPrms.tau = fluorTau;

% fit parameters for decay - 
fts = fittype('peak*exp(-x/tau)');
fos = fitoptions(fts);
fos.StartPoint = [1, 0.05];
fos.Lower = [0 0];
fos.Upper = [2 1000];
fitPrms.fitType = fts;
fitPrms.fitOptions = fos;
fitPrms.startFitOffset = 0.005; % 5ms

% optimizationParameters - Range
restRange = [1e-9 500e-9];

betaRange = [500 3000];
iCalciumAmpRange = [0 1e-6];
konRange = [1e5 5e11];
kdRange = [10e-9 1];
btRange = [100e-9 1];
% To go in optimization engine
iPrm = [50e-9, 1400, 7e-12, 5e8, 100e-6, 1e-3];
lowerBound = [restRange(1), betaRange(1), iCalciumAmpRange(1), konRange(1), kdRange(1), btRange(1)];
upperBound = [restRange(2), betaRange(2), iCalciumAmpRange(2), konRange(2), kdRange(2), btRange(2)];
% Optimization Parameters
lsqOptions = optimoptions('lsqnonlin','OptimalityTolerance',1e-20,'FunctionTolerance',1e-20,'StepTolerance',1e-10,...
    'MaxFunctionEvaluations',1e3,'Display','iter-detailed');


% Do optimization
objective = @(optPrms) optimizeSystem(optPrms,systemPrms,fluorPrms,odePrms,fitPrms,odeOptions);
sol = lsqnonlin(objective, iPrm, lowerBound, upperBound, lsqOptions);


%% ODEs used above

function cost = optimizeSystem(optPrms,systemPrms,fluorPrms,odePrms,fitPrms,odeOptions)
    NO = sum(fluorPrms.useBufferIdx); % How many objectives are we running
    objIdx = find(fluorPrms.useBufferIdx);
    
    referenceValues = cat(1, fluorPrms.peak(fluorPrms.useBufferIdx)', fluorPrms.tau(fluorPrms.useBufferIdx)');
    estimateValues = zeros(NO,2);
    
    for no = 1:NO
        cObjIdx = objIdx(no);
        
        currentOptPrms.rest= optPrms(1);
        currentOptPrms.beta = optPrms(2);
        currentOptPrms.iCalciumAmp = optPrms(3);
        currentOptPrms.kon = optPrms(4);
        currentOptPrms.kd = optPrms(5);
        currentOptPrms.bt = optPrms(6);
        
        currentFluorPrms.kon = fluorPrms.kon; % For now this is a single value
        currentFluorPrms.kd = fluorPrms.kd(cObjIdx);
        currentFluorPrms.bt = fluorPrms.bt(cObjIdx);
        
        cKappaEB = optPrms(6) / (optPrms(1) + optPrms(5));
        cKappaFB = currentFluorPrms.bt / (optPrms(1) + currentFluorPrms.kd);
        cu0 = [optPrms(1) optPrms(1)*cKappaEB optPrms(1)*cKappaFB];
        
        odeFunction = @(t,y) endogenousAndFluorophores(t,y,systemPrms,currentFluorPrms,currentOptPrms);
        [cTime, cVals] = ode23s(odeFunction, odePrms.tspan, cu0, odeOptions);
        estimateValues(no,1) = max(cVals(:,1)); % Estimate of peak calcium 
        [~,fluorPkIdx] = max(cVals(:,3));
        fitIdx = find(cTime>= (cTime(fluorPkIdx)+fitPrms.startFitOffset),1,'first');
        expTime = cTime(fitIdx:end)-cTime(fitIdx);
        expValues = cVals(fitIdx:end,3)-cVals(end,3);
        currentFit = fit(expTime,expValues/expValues(1),fitPrms.fitType,fitPrms.fitOptions);
        estimateValues(no,2) = 1000*currentFit.tau;
    end
    cost = (estimateValues(:) - referenceValues) ./ referenceValues;
end

% ODE for endogenous + fluorescent buffer and optimization of system prms
function dydt = endogenousAndFluorophores(t,y,systemPrms,fluorPrms,optPrms)
    % Calcium current stuff
    restCurrent = optPrms.beta * optPrms.rest;
    curr = restCurrent + ica(t,optPrms.iCalciumAmp)/(2*96485*systemPrms.v);
    extrusion = optPrms.beta * y(1);

    % Endogenous Buffer Reaction 
    endogAssociation = y(1)*(optPrms.bt - y(2))*optPrms.kon;
    endogDissociation = y(2)*optPrms.kon*optPrms.kd;
    endogExchange = endogAssociation - endogDissociation;

    % Fluorescent Buffer Reaction
    fluorAssociation = y(1)*(fluorPrms.bt - y(3))*fluorPrms.kon;
    fluorDissociation = y(3)*fluorPrms.kon*fluorPrms.kd;
    fluorExchange = fluorAssociation - fluorDissociation;

    %Output
    dydt(1,1) = curr - extrusion - endogExchange - fluorExchange;
    dydt(2,1) = endogExchange;
    dydt(3,1) = fluorExchange;
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












































































