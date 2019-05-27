function out = simulateBarrage(timeprm,synprm,stimprm,recprm,anprm)


% Keep these in output
out.timeprm = timeprm;
out.synprm = synprm;
out.stimprm = stimprm;
out.recprm = recprm;
out.anprm = anprm;

% For making the conductances
alpha = @(t, rise, fall) (t>=0).*(exp(-t/fall) - exp(-t/rise)); % only positive t's...


%% -- get parameters --

% Time
dt = timeprm.dt; % ms
T = timeprm.T; % ms
out.tvec = 0:dt:T; 
NT = length(out.tvec); % number of samples

% Synaptic Parameters
excRise = synprm.excRise; % ms
excFall = synprm.excFall; % ms
excAmp = synprm.excAmp; % siemens
excRev = synprm.excRev; % volts
excFreq = synprm.excFreq; % 1/s
excNumber = excFreq * T/1000; % number of inputs

inhRise = synprm.inhRise; % ms
inhFall = synprm.inhFall; % ms
inhAmp = synprm.inhAmp; % siemens
inhRev = synprm.inhRev; % volts
inhFreq = synprm.inhFreq; % 1/sec
inhNumber = inhFreq * T/1000; % Number of inputs

% Stimulation Parameters
holdVoltage = stimprm.holdVoltage; % volts
modulationDepth = stimprm.modulationDepth; % volts
modulationPeriod = stimprm.modulationPeriod; % ms

% Recording Parameters
noiseAmplitude = recprm.noiseAmplitude; % amperes (standard deviation)

% Analysis Parameters
numCycles = anprm.numCycles; % number of modulation cycles to average over


%% -- generate conductances --

out.eTrain = randi(T,excNumber, 1); % time of each conductance
out.iTrain = randi(T,inhNumber, 1);
out.eConductance = zeros(1,NT); % excitatory conductance
out.iConductance = zeros(1,NT); % inhibitory conductance
for ne = 1:excNumber
    currCond = excAmp * alpha(out.tvec-eTrain(ne),excRise,excFall);
    currCond(isnan(currCond))=0;
    out.eConductance = eConductance + currCond;
end
for ni = 1:inhNumber
    currCond = inhAmp * alpha(out.tvec-iTrain(ni),inhRise,inhFall);
    currCond(isnan(currCond)) = 0;
    out.iConductance = iConductance + currCond;
end


%% -- recording simulation --

out.voltage = holdVoltage + modulationDepth/2*sin(2*pi*out.tvec/modulationPeriod); % voltage command
out.eCurrent = eConductance .* (voltage - excRev); % Excitatory Current
out.iCurrent = iConductance .* (voltage - inhRev); % Inhibitory Current
out.synCurrent = eCurrent + iCurrent; % Total Synaptic Current
out.nCurrent = noiseAmplitude*randn(1,NT); % Noise current
out.totalCurrent = synCurrent + nCurrent; % Total Measured Current


%% -- perform analysis --

aWindowTime = numCycles * modulationPeriod; % Time of analysis window
aWindowSamples = aWindowTime/dt; % Samples of analysis window
aWindowStart = 1:aWindowSamples:NT; % Start of each window in samples
out.aWindowCenter = mean([out.tvec(aWindowStart(1:end-1));out.tvec(aWindowStart(2:end))],1); % Center of each window in time
NAW = length(aWindowStart); % Number of analysis windows

out.aResult = zeros(aWindowSamples,NAW-1,2); 
out.aLine = zeros(2,NAW-1);
out.estConductance = zeros(2,NAW-1); % Estimate of conductance [gExc; gInh]
out.cycleConductance = zeros(2,NAW-1); % Average of true conductance per cycle
out.residual = zeros(2,NAW-1); % Residual (based on average of each analysis window)
for naw = 1:NAW-1
    cSamples = aWindowStart(naw):aWindowStart(naw+1)-1;
    
    % Results --
    out.aResult(:,naw,1) = out.totalCurrent(cSamples);
    out.aResult(:,naw,2) = out.voltage(cSamples);
    out.aLine(:,naw) = [out.aResult(:,naw,2) ones(aWindowSamples,1)] \ out.aResult(:,naw,1);
    
    % Convert to Conductances
    cgSyn = aLine(1,naw);
    ceSyn = -aLine(2,naw) / cgSyn;
    cgi = (cgSyn*ceSyn - cgSyn*excRev) / (inhRev - excRev);
    cge = cgSyn - cgi;
    out.estConductance(1,naw) = cge;
    out.estConductance(2,naw) = cgi;
    
    % Compute Residual
    out.cycleConductance(:,naw) = [mean(eConductance(cSamples)); mean(iConductance(cSamples))];
    out.residual(:,naw) = out.estConductance(:,naw) - out.cycleConductance(:,naw);
end

%% -- generate summary statistics --

out.rmsError = sqrt(mean(residual.^2),2); % rms error
[~,excGOF] = fit(cycleConductance(1,:)',out.estConductance(1,:)','poly1'); % linear fit to get R^2
[~,inhGOF] = fit(cycleConductance(2,:)',out.estConductance(2,:)','poly1'); % linear fit to get R^2
out.excr2 = excGOF.rsquare;
out.inhr2 = inhGOF.rsquare;












    









