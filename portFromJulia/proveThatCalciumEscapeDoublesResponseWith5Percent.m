


hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/portFromJulia';
fpath = '/Users/landauland/Documents/Julia/calciumBuffering/latexDiscussion/images';

% Spine Parameters
systemPrms.beta=1400.0;
systemPrms.rest=5.0e-8;
systemPrms.amp=3.3e-12;
systemPrms.ks=20.0;
systemPrms.v=(4/3)*pi*0.7^3*1e-15;

u0 = systemPrms.rest;
tspan = [0.001,0.01];
dt = 0.000001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',dt);

NA = 15;
amp = logspace(-12,-10.5,NA);

instantSols = zeros(NT,NA);
escapeSols = zeros(NT,NA);


koff = 1.06e12;
kon = koff/175e-6;

for na = 1:NA
    systemPrms.amp = amp(na);
    isol = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms),tspan,u0,odeOptions);
    instantSols(:,na) = deval(isol,tvec);
    
    btEstimate = max(instantSols(:,na))*15.5*2;
    prms = [kon koff btEstimate];
    kappa = btEstimate / (systemPrms.rest + koff/kon);
    fprintf('%d \n',kappa);
    e0 = [systemPrms.rest systemPrms.rest*kappa];
    esol = ode23s(@(t,y) massActionDynamics(t,y,systemPrms,prms),tspan,e0,odeOptions);
    cvals = deval(esol,tvec);
    escapeSols(:,na) = cvals(1,:);
end

f = figure(1); clf;
subplot(1,3,1);
plotvc(tvec,instantSols,f,{'''linewidth''','2'});
subplot(1,3,2);
plotvc(tvec, escapeSols,f,{'''linewidth''','2'});
subplot(1,3,3);
plotvc(tvec,escapeSols./instantSols,f,{'''linewidth''','2'});


%% Functions
function out = ica(t,amp)
    % function giving calcium current with fixed temporal properties and
    % user-defined amplitude
    current = amp*exp(-(t-2e-3).^2 / (0.55e-3)^2);
    offset = amp*exp(-(1e-3-2e-3).^2 / (0.55e-3)^2);
    adjustedCurrent = (current - offset) .* (t>=1e-3).*(t<=3e-3);
    out = adjustedCurrent .* (adjustedCurrent>=0);
end

% Gives solution with instant buffer
function dy = instantBufferDynamics(t,y,prm)
    % Current/Extrusion Term
    restCurrent = prm.beta/prm.ks * prm.rest;
    curr = restCurrent + ica(t,prm.amp)/(96485*prm.v*prm.ks);
    extrusion = prm.beta/prm.ks * y(1);
    % Derivative
    dy = curr - extrusion;
end

% Gives solution with mass-action buffer
function dy = massActionDynamics(t,y,systemPrms,prms)
    % Parameters: [v,ks,beta,rest,amp]
    % Kon,Koff,Bt are optimization variables
    
    % y(1) = cafree
    % y(2) = boundBuffer
    % y(3) = secondBoundBuffer
    
    % Current/Extrusion Term
    restCurrent = systemPrms.beta * systemPrms.rest;
    curr = restCurrent + ica(t,systemPrms.amp)/(96485*systemPrms.v);
    extrusion = systemPrms.beta * y(1);
    
    % Association/Dissociation with quasi-instant buffer
    association = y(1)*(prms(3)-y(2))*prms(1);
    dissociation = y(2)*prms(2);
    bufferExchange = association - dissociation;
    
    % Derivative
    dy(1,1) = curr - extrusion - bufferExchange;
    dy(2,1) = bufferExchange;
end