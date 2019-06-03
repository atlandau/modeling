function out = vcStimulation(tprm,exc,inh,stim,cellprm)
% adds a line to convert holdVoltage to cell voltage that is phase shifted
% and damped by the cell's membrane

% Log Normal Inline
lnprm = @(mn,vr) [log((mn^2)/sqrt(vr+mn^2)), sqrt(log(vr/(mn^2)+1))];
alpha = @(t, rise, fall) double(t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's...

% This flips the polarity of the sine wave modulation each run.
% It's gated by stim.modulationShift.
persistent phaseSwitch
if isempty(phaseSwitch)
    phaseSwitch = true;
end
phaseSwitch = ~phaseSwitch;

% Time Parameters
out.tvec = 0:tprm.dt:tprm.T; %vector (ms)
NT = length(out.tvec);


% Synaptic Parameters
enPrms = lnprm(exc.numberMean,exc.numberVar); % Excitatory Number - LN Parameters
edPrms = lnprm(exc.delayMean,exc.delayVar); % Excitatory Delay - LN Parameters
inPrms = lnprm(inh.numberMean,inh.numberVar); % Inhibitory Number - LN Parameters
idPrms = lnprm(inh.delayMean,inh.delayVar); % Inhibitory Delay - LN Parameters

% Generate Conductances
eNumber = round(lognrnd(enPrms(1),enPrms(2)));
iNumber = round(lognrnd(inPrms(1),inPrms(2)));
% Generate timing trains, resolve to dt, add hard delay and baselineWindow
eTrain = tprm.T/2 + exc.hardDelay + round(lognrnd(edPrms(1),edPrms(2),[eNumber 1])/tprm.dt)*tprm.dt;
iTrain = tprm.T/2 + inh.hardDelay + round(lognrnd(idPrms(1),idPrms(2),[iNumber 1])/tprm.dt)*tprm.dt;
excConds = zeros(eNumber,NT);
inhConds = zeros(iNumber,NT);
for ne = 1:eNumber
    % Generate each conductance train and put in array 
    cConductance = exc.excAmp * alpha(out.tvec-eTrain(ne),exc.rise,exc.fall);
    cConductance(isnan(cConductance))=0; % Values too large become nans, make them 0
    excConds(ne,:) = cConductance;
end
for ni = 1:iNumber
    cConductance = inh.inhAmp * alpha(out.tvec-iTrain(ni),inh.rise,inh.fall);
    cConductance(isnan(cConductance)) = 0; % Values too large become nans, make them 0
    inhConds(ni,:) = cConductance;
end
out.eConductance = sum(excConds,1); % Sum up all conductances from each synapse
out.iConductance = sum(inhConds,1);

% Create synaptic parameters for simulation
syn.ee = exc.rev;
syn.ei = inh.rev;
syn.tvec = out.tvec;
syn.ge = out.eConductance;
syn.gi = out.iConductance;
syn.openCircuit = cellprm.openCircuit;


% Stimulation
varyPhase = pi*stim.modulationShift*phaseSwitch + rand*2*pi*(stim.modulationShift==-1);
cellprm.vc = @(t) stim.vHold + stim.modulationDepth/2 * sin(2*pi*t/stim.modulationPeriod + varyPhase); % Sine wave mod
iState = cellprm.vc(0);
out.holdVoltage = cellprm.vc(out.tvec);
[tvecCheck,out.cellVoltage] = ...
    eulerapp(@(t,v) vcSynapticDiffEQ(t,v,cellprm,syn),[0 tprm.T],iState,tprm.dt,1,0); % Generate true cell voltage

if ~isequal(tvecCheck,out.tvec)
    error('Time vector generated incorrectly in eulerapp...');
end
out.cellVoltage = out.cellVoltage(:)';

% Estimate Voltage with circuit analysis
holdingCenter = (cellprm.em*cellprm.rs + stim.vHold*cellprm.rm) / (cellprm.rs + cellprm.rm);
omega = 2*pi/stim.modulationPeriod;
realImpedance = (cellprm.rm + cellprm.rs)/cellprm.rm;
imagImpedance = omega * cellprm.rs * cellprm.cm;
magResponse = stim.modulationDepth / sqrt(realImpedance^2 + imagImpedance^2);
timeDelay = atan(imagImpedance/realImpedance)/(2*pi) * stim.modulationPeriod;
out.estimateVoltage = holdingCenter + ...  
    magResponse/2*(2/stim.modulationDepth*(cellprm.vc(out.tvec-timeDelay)-stim.vHold));


% Compute Currents
out.eCurrent = syn.ge .* (out.cellVoltage - exc.rev); % Excitatory Current
out.iCurrent = syn.gi .* (out.cellVoltage - inh.rev); % Inhibitory Current
out.rCurrent = (out.cellVoltage - cellprm.em) / cellprm.rm; % resistive current
out.cCurrent = cellprm.cm * diff(out.cellVoltage)/tprm.dt; % capacitive current
out.cCurrent(end+1) = out.cCurrent(end); % add in relatively accurate last value
out.synCurrent = out.eCurrent + out.iCurrent; % Total Synaptic Current
out.nCurrent = stim.noiseAmplitude*randn(1,length(out.tvec)); % Noise current
out.totalCurrent = out.synCurrent + out.nCurrent + out.rCurrent + out.cCurrent; % Total Measured Current

% Do subtraction
cyclesSkipped = 3;
samplesPerCycle = stim.modulationPeriod/tprm.dt;
baseCurrentNumberCycles = floor(tprm.T/2 / stim.modulationPeriod);
baseCurrent = out.totalCurrent(cyclesSkipped*samplesPerCycle+1:baseCurrentNumberCycles*samplesPerCycle);
out.averageBaseCurrent = mean(reshape(baseCurrent,samplesPerCycle,baseCurrentNumberCycles-cyclesSkipped),2)';
numberFullCycles = floor(NT/samplesPerCycle);
remainingSamples = NT - numberFullCycles*samplesPerCycle;
out.subtractCurrent = [repmat(out.averageBaseCurrent,1,numberFullCycles), out.averageBaseCurrent(1:remainingSamples)];
out.estimateCurrent = out.totalCurrent - out.subtractCurrent;

% Analysis Cycles
aWindowTime = stim.aCycles * stim.modulationPeriod; % Time of analysis window
aWindowSamples = aWindowTime/tprm.dt; % Samples of analysis window

if isnumeric(stim.method)
    % Lazy but this works...
    downSample = round(aWindowTime*stim.method/tprm.dt);
    NAW = NT-(aWindowSamples-1);
    aWindowStart = 1:NAW;
    aWindowStart = aWindowStart(1:downSample:end);
    NAW = length(aWindowStart);
else
    switch stim.method
        case 'simple'
            NAW = floor(NT/aWindowSamples);
            aWindowStart = 1:aWindowSamples:aWindowSamples*NAW;
        case 'interleaved'
            NAW = NT-(aWindowSamples-1);
            aWindowStart = 1:NAW; % Start of each window in samples
    end
end

out.aWindowCenter = out.tvec(aWindowStart)+aWindowTime/2;
out.aResult = zeros(aWindowSamples,NAW,2); 
out.aLine = zeros(2,NAW);
out.estConductance = zeros(2,NAW); % Estimate of conductance [gExc; gInh]
out.cycleConductance = zeros(2,NAW); % Average of conductance per analysis cycle
out.residual = zeros(2,NAW); % Residual (based on average of each analysis window)
for naw = 1:NAW    
    cSamples = aWindowStart(naw):aWindowStart(naw)+aWindowSamples-1;
    
    % Results --
    out.aResult(:,naw,1) = out.estimateCurrent(cSamples);
    out.aResult(:,naw,2) = out.estimateVoltage(cSamples);
    out.aLine(:,naw) = [out.aResult(:,naw,2) ones(aWindowSamples,1)] \ out.aResult(:,naw,1); % Linear Regression
    
    % Convert to Conductances (Wehr & Zador, 2003, Nature)
    % Note that this allows negative conductance... could fix?
    cgSyn = out.aLine(1,naw);
    ceSyn = -out.aLine(2,naw) / cgSyn;
    cgi = (cgSyn*ceSyn - cgSyn*exc.rev) / (inh.rev - exc.rev);
    cge = cgSyn - cgi;
    out.estConductance(1,naw) = cge;
    out.estConductance(2,naw) = cgi;
    
    % Compute Residual
    out.cycleConductance(:,naw) = [mean(out.eConductance(cSamples)); mean(out.iConductance(cSamples))];
    out.residual(:,naw) = out.estConductance(:,naw) - out.cycleConductance(:,naw);
end

out.rmsError = sqrt(mean((1e9*out.residual(round(end/2):end)).^2,2)); % rms error
[~,excGOF] = fit(1e9*out.cycleConductance(1,round(end/2):end)',1e9*out.estConductance(1,round(end/2):end)','poly1'); % linear fit to get R^2
[~,inhGOF] = fit(1e9*out.cycleConductance(2,round(end/2):end)',1e9*out.estConductance(2,round(end/2):end)','poly1'); % linear fit to get R^2
out.excr2 = excGOF.rsquare;
out.inhr2 = inhGOF.rsquare;








    







