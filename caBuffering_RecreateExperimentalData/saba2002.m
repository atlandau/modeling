

hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/caBuffering_RecreateExperimentalData';
fpath = fullfile(hpath,'figures');

% Spine Parameters
systemPrms.beta=1400.0;
systemPrms.rest=5.0e-8;
systemPrms.amp=6e-12;
systemPrms.ks=20.0;
systemPrms.v=(4/3)*pi*0.7^3*1e-15;

u0 = systemPrms.rest;
tspan = [0,0.1];
dt = 0.000001;
tvec = tspan(1):dt:tspan(2);

% Instant Solution
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',dt);
[trueTime,trueSol] = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms),tspan,u0,odeOptions);
plot(trueTime,trueSol);


%% Three Different Models
% model types:
bufferType = {'lipid','atp','cam'};








%% ODEs used above


% ODE for instant, never-saturated buffer
function dy = instantBufferDynamics(t,y,prm)
    % Current/Extrusion Term
    restCurrent = prm.beta/prm.ks * prm.rest;
    curr = restCurrent + ica(t,prm.amp)/(2*96485*prm.v*prm.ks);
    extrusion = prm.beta/prm.ks * y(1);
    % Derivative
    dy = curr - extrusion;
end

% ODE for endogenous buffer only
function dy = massActionDynamics(t,y,systemPrms,prms)
    % SystemParameters: [v,ks,beta,rest,amp]
    % prms: [Kon,Koff,Bt] of endogenous buffer
    % y(1) = cafree
    % y(2) = boundBuffer
    
    % Current/Extrusion Term
    restCurrent = systemPrms.beta * systemPrms.rest;
    curr = restCurrent + ica(t,systemPrms.amp)/(2*96485*systemPrms.v);
    extrusion = systemPrms.beta * y(1);
    % Association/Dissociation with quasi-instant buffer
    association = y(1)*(prms(3)-y(2))*prms(1);
    dissociation = y(2)*prms(2);
    bufferExchange = association - dissociation;
    % Derivative
    dy(1,1) = curr - extrusion - bufferExchange;
    dy(2,1) = bufferExchange;
end

% ODE for endogenous buffer and fluorescent buffer
function dy = twoBufferDynamics(t,y,systemPrms,prms)
    % SystemParameters: [v,ks,beta,rest,amp]
    % prms: [eKon,eKoff,eBt,fKon,fKoff,fBt] e/f::endogenous/fluorescent
    % y(1) = cafree
    % y(2) = boundEndogenousBuffer
    % y(3) = boundFluorescentBuffer
    
    % Current/Extrusion Term
    restCurrent = systemPrms.beta * systemPrms.rest;
    curr = restCurrent + ica(t,systemPrms.amp)/(2*96485*systemPrms.v);
    extrusion = systemPrms.beta * y(1);
    % Association/Dissociation with endogenous buffer
    eAssociation = y(1)*(prms(3)-y(2))*prms(1);
    eDissociation = y(2)*prms(2);
    eExchange = eAssociation - eDissociation;
    % Association/Dissociation with fluorescent buffer
    fAssociation = y(1)*(prms(6)-y(3))*prms(4);
    fDissociation = y(3)*prms(5);
    fExchange = fAssociation - fDissociation;
    % Derivative
    dy(1,1) = curr - extrusion - eExchange - fExchange;
    dy(2,1) = eExchange;
    dy(3,1) = fExchange;
end








