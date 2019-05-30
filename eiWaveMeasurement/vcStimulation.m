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
enPrms = lnprm(exc.excNumberMean,exc.excNumberVar); % Excitatory Number - LN Parameters
edPrms = lnprm(exc.excDelayMean,exc.excDelayVar); % Excitatory Delay - LN Parameters
inPrms = lnprm(inh.inhNumberMean,inh.inhNumberVar); % Inhibitory Number - LN Parameters
idPrms = lnprm(inh.inhDelayMean,inh.inhDelayVar); % Inhibitory Delay - LN Parameters


% Generate Conductances
eNumber = round(lognrnd(enPrms(1),enPrms(2)));
iNumber = round(lognrnd(inPrms(1),inPrms(2)));
% Generate timing trains, resolve to dt, add hard delay
eTrain = exc.excHardDelay + round(lognrnd(edPrms(1),edPrms(2),[eNumber 1])/tprm.dt)*tprm.dt;
iTrain = inh.inhHardDelay + round(lognrnd(idPrms(1),idPrms(2),[iNumber 1])/tprm.dt)*tprm.dt;
excConds = zeros(eNumber,NT);
inhConds = zeros(iNumber,NT);
for ne = 1:eNumber
    % Generate each conductance train and put in array 
    cConductance = exc.excAmp * alpha(out.tvec-eTrain(ne),exc.excRise,exc.excFall);
    cConductance(isnan(cConductance))=0; % Values too large become nans, make them 0
    excConds(ne,:) = cConductance;
end
for ni = 1:iNumber
    cConductance = inh.inhAmp * alpha(out.tvec-iTrain(ni),inh.inhRise,inh.inhFall);
    cConductance(isnan(cConductance)) = 0; % Values too large become nans, make them 0
    inhConds(ni,:) = cConductance;
end
out.eConductance = sum(excConds,1); % Sum up all conductances from each synapse
out.iConductance = sum(inhConds,1);

% Stimulation
cellprm.vc = @(t) stim.vHold + ... % adjust to hold voltage
    stim.modulationDepth/2 * ... % set height 
    sin(2*pi*t/stim.modulationPeriod + phaseSwitch*pi*stim.modulationShift); % create sine wave
cellprm.em = cellprm.vc(0);
iState = cellprm.em;
out.holdVoltage = cellprm.vc(out.tvec);
[tvecCheck,out.cellVoltage] = eulerapp(@(t,v) vcdiffeq(t,v,cellprm),[0 tprm.T],iState,tprm.dt,1,0); % Generate true cell voltage
if ~isequal(tvecCheck,out.tvec)
    error('Time vector generated incorrectly in eulerapp...');
end
out.cellVoltage = out.cellVoltage(:)';

out.eCurrent = out.eConductance .* (out.cellVoltage - exc.excRev); % Excitatory Current
out.iCurrent = out.iConductance .* (out.cellVoltage - inh.inhRev); % Inhibitory Current
synCurrent = out.eCurrent + out.iCurrent; % Total Synaptic Current
nCurrent = stim.noiseAmplitude*randn(1,length(out.tvec)); % Noise current
out.totalCurrent = synCurrent + nCurrent; % Total Measured Current

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
aResult = zeros(aWindowSamples,NAW,2); 
out.aLine = zeros(2,NAW);
out.estConductance = zeros(2,NAW); % Estimate of conductance [gExc; gInh]
out.cycleConductance = zeros(2,NAW); % Average of conductance per analysis cycle
out.residual = zeros(2,NAW); % Residual (based on average of each analysis window)
for naw = 1:NAW    
    cSamples = aWindowStart(naw):aWindowStart(naw)+aWindowSamples-1;
    
    % Results --
    aResult(:,naw,1) = out.totalCurrent(cSamples);
    aResult(:,naw,2) = out.holdVoltage(cSamples);
    out.aLine(:,naw) = [aResult(:,naw,2) ones(aWindowSamples,1)] \ aResult(:,naw,1); % Linear Regression
    
    % Convert to Conductances (Wehr & Zador, 2003, Nature)
    % Note that this allows negative conductance... could fix?
    cgSyn = out.aLine(1,naw);
    ceSyn = -out.aLine(2,naw) / cgSyn;
    cgi = (cgSyn*ceSyn - cgSyn*exc.excRev) / (inh.inhRev - exc.excRev);
    cge = cgSyn - cgi;
    out.estConductance(1,naw) = cge;
    out.estConductance(2,naw) = cgi;
    
    % Compute Residual
    out.cycleConductance(:,naw) = [mean(out.eConductance(cSamples)); mean(out.iConductance(cSamples))];
    out.residual(:,naw) = out.estConductance(:,naw) - out.cycleConductance(:,naw);
end

out.rmsError = sqrt(mean((1e9*out.residual).^2,2)); % rms error
[~,excGOF] = fit(1e9*out.cycleConductance(1,:)',1e9*out.estConductance(1,:)','poly1'); % linear fit to get R^2
[~,inhGOF] = fit(1e9*out.cycleConductance(2,:)',1e9*out.estConductance(2,:)','poly1'); % linear fit to get R^2
out.excr2 = excGOF.rsquare;
out.inhr2 = inhGOF.rsquare;








    







