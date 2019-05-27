function gridSearch_barrage(modDepthIDX,modPeriodIDX,numCyclesIDX)


% do outer loop on modulationDepth and modulationPeriod
% everything else goes into inner loop

%% Grid Parameters
excFall = [1 3 5 10 15]; % ms
excAmp = 1e-9*[0.1 0.5 1 3 5 10]; % siemens
excFreq = [10 25 50 100 200]; % 1/sec

inhFall = [3 5 10 15 30]; % ms
inhAmp = 1e-9*[0.1 0.5 1 3 5 10]; % siemens
inhFreq = [10 25 50 100 200]; % 1/sec

modulationDepth = 1e-3*[0.1 0.25 0.5 1 2.5 5 10 25 50]; % volts
modulationPeriod = [0.1 0.25 0.5 1 2.5 5 10]; % ms

noiseAmplitude =  1e-12*[1 2.5 5 10 25]; % amperes

numCycles = 1:10; % cycles of modulation to average over


%% Inner Loop Handling
NEF = length(excFall);
NEA = length(excAmp);
NEH = length(excFreq);
NIF = length(inhFall);
NIA = length(inhAmp);
NIH = length(inhFreq);
NNA = length(noiseAmplitude);


%% Create Parameter Structures

for nef = 1:NEF
    for nea = 1:NEA
        for neh = 1:NEH
            for nif = 1:NIF
                for nia = 1:NIA
                    for nih = 1:NIH
                        for nna = 1:NNA
                            % Time Parameters
                            timeprm.dt = 0.01; %ms
                            timeprm.T = 10000; % ms

                            % Synaptic Parameters
                            synprm.excRise = 0.3; % ms
                            synprm.excFall = excFall(neh); 
                            synprm.excAmp = excAmp(nea);
                            synprm.excRev = 0e-3; %volts
                            synprm.excFreq = excFreq(neh);

                            synprm.inhRise = 2.5; % ms
                            synprm.inhFall = inhFall(nif);
                            synprm.inhAmp = inhAmp(nia);
                            synprm.inhRev = -70e-3; % volts
                            synprm.inhFreq = inhFreq(nih);
                            
                            % Stimulation Parameters
                            stimprm.holdVoltage = -35e-3; % volts
                            stimprm.modulationDepth = modulationDepth(modDepthIDX); % outer loop 
                            stimprm.modulationPeriod = modulationPeriod(modPeriodIDX); % outer loop
                            
                            % Recording Parameters (include capacitance)
                            recprm.noiseAmplitude = noiseAmplitude(nna); 
                            
                            % Analysis Parameteres
                            anprm.numCycles = numCycles(numCyclesIDX); % outer loop
                            
                            % Do simulation
                            out = simulateBarrage(timeprm,synprm,stimprm,recprm,anprm);
                            
                            name = sprintf('Barrage_md%d_mp%d_nc%d_ef%d_ea%d_eh%d_if%d_ia%d_ih%d_na%d.mat',...
                                modDepthIDX,modPeriodIDX,numCyclesIDX,...
                                nef,nea,neh,nif,nia,nih,nna);
                            savepath = fullfile('~/synapticBarrage_190527',name);
                            save(savepath,'out');
                        end
                    end
                end
            end
        end
    end
end






















    


