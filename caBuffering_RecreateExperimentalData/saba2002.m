
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
tspan = [0,0.1];
dt = 0.000001;
tvec = tspan(1):dt:tspan(2);

% Instant Solution
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',dt);
[trueTime,trueSol] = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms,iCalciumAmp),tspan,u0,odeOptions);
plot(trueTime,trueSol);


%% Three Different Models
% model types:
bufferType = {'lipid','atp','cam'};

% fluorescent buffers
names = {'fluo5f','ogb1','fluo4'};
useBufferIdx = [0 1 1];
% fluorescent parameters
fluorKd = [800e-9, 205e-9, 300e-9];
fluorKon = [5e8, 5e8, 5e8];
fluorConcentrations = {300e-6, [20e-6, 50e-6, 100e-6], 20e-6};
fluorPeak = {77e-9, [400e-9, 162.7e-9, 68.9e-9], 474.5e-9};
fluorTau = {100, [30.7 67.5 141.1], 20.77}; 
% fluorescent parameter structure
fluorPrms.idx = useBufferIdx;
fluorPrms.kd = fluorKd;
fluorPrms.kon = fluorKon;
fluorPrms.concentration = fluorConcentrations;
fluorPrms.peak = fluorPeak;
fluorPrms.tau = fluorTau;

% Calcium Amplitudes
iCalciumAmplitudes = [7, 21, 70] * 1e-12; 


% Now - endogenous buffer prms
% prm array that goes into the objective function
% -- is going to have systemPrms.beta in it
% -- and will have for each endogenous buffer type various different prms


% CaM - dydt = camODE(t,y,iCalciumAmp,systemPrms,camPrms,fluorPrms)









%% - - - example lsqnonlin

konStart = 5e9;
kdStart = 200e-6;
iPrm = [konStart; kdStart];
lowerBound = [konRange(1); kdRange(1)];
upperBound = [konRange(2); kdRange(2)];

bestOptions = zeros(2,NBT);
bestTrajectories = zeros(NT,NA,NBT);

lsqOptions = optimoptions('lsqnonlin','OptimalityTolerance',1e-20,'FunctionTolerance',1e-20,'StepTolerance',1e-10);
for nbt = 1:NBT
    fprintf('%d/%d...\n',nbt,NBT);
    
    objective = @(prms) objFcn(prms,systemPrms,odeOptions,tvec,instantSols,amp,bt(nbt),tspan);
    sol = lsqnonlin(objective, iPrm, lowerBound, upperBound, lsqOptions);
    bestOptions(:,nbt) = sol;
    
    for na = 1:NA
        systemPrms.amp = amp(na);
        kappa = bt(nbt) / (systemPrms.rest + bestOptions(2));
        u0 = [systemPrms.rest systemPrms.rest*kappa];
        bsol = ode23s(@(t,y) massActionDynamics_FixBT(t,y,systemPrms,bt(nbt),bestOptions(:,na)), tspan,u0,odeOptions);
        bvals = deval(bsol,tvec);
        bestTrajectories(:,na,nbt) = bvals(1,:);
    end
end


%% ODEs used above

function cost = costCaMFluors(prms,fluorPrms,systemPrms,odeOptions,tvec,iCalciumAmplitudes,tSpan)
    
end



function cost = objFcn(prms,systemPrms,odeOptions,tvec,instantSols,amplitudes,bt,tSpan)
    testSols = zeros(length(tvec),length(amplitudes));
    for na = 1:length(amplitudes)
        systemPrms.amp = amplitudes(na);
        kappa = bt / (systemPrms.rest + prms(2)/prms(1));
        u0 = [systemPrms.rest systemPrms.rest*kappa];
        curSol = ode23s(@(t,y) massActionDynamics_FixBT(t,y,systemPrms,bt,prms), tSpan,u0,odeOptions);
        curVals = deval(curSol,tvec);
        testSols(:,na) = curVals(1,:);
    end
    cost = testSols(:) - instantSols(:);
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





function cost = costLipidFluors(prms,fluorPrms,systemPrms,odeOptions,tvec,iCalciumAmplitudes,tSpan)
end

function cost = costAtpFluors(prms,fluorPrms,systemPrms,odeOptions,tvec,iCalciumAmplitudes,tSpan)
end







































































