
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

%% -- vary Kon and Concentration of Endogenous

tspan = [0,0.02];
odeOptions = odeset('RelTol',1e-8,'AbsTol',1e-8);%,'MaxStep',0.02);
dt = 0.0001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

spineVolume = 1e-15; % 1fL
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
fluo4ff.totalConcentration = 100e-6;
fluo4ff.kd = 9e-6;

fluorProperties = fluo4ff;
fluorKappa = kappaFunc(fluorProperties,restConcentration);

% Endogenous Buffer Properties
numOnRate = 10;
numConcentration = 20;
onRateRange = [4 12];
concentrationRange = [-7 -1];
onRate = logspace(onRateRange(1),onRateRange(2),numOnRate);
concentration = logspace(concentrationRange(1),concentrationRange(2),numConcentration);
kd = concentration/20; % keep Kd optimal with concentration

msg = '';
values = zeros(NT,3,numOnRate,numConcentration);
for non = 1:numOnRate
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf('On Rate %d/%d ... \n',non,numOnRate);
    fprintf(1,msg);
    bufferProperties.onRate = onRate(non);
    for nc = 1:numConcentration
        bufferProperties.totalConcentration = concentration(nc);
        bufferProperties.kd = kd(nc);
        iState = restConcentration * [1 kappaFunc(bufferProperties,restConcentration) fluorKappa];
            
        prmODE = @(t,y) endogenousAndFluorescent(t,y,sysProperties,bufferProperties,fluorProperties);
        csol = ode23s(prmODE,tspan,iState,odeOptions);
        values(:,:,non,nc) = deval(csol,tvec)';
    end
end






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



































































