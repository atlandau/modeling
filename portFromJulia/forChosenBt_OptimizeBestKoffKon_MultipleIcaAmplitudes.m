

hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/portFromJulia';
fpath = '/Users/landauland/Documents/Julia/calciumBuffering/latexDiscussion/images';

amp = [3.3e-12 10e-12 35e-12];
NA = length(amp);

% Spine Parameters
systemPrms.beta=1400.0;
systemPrms.rest=5.0e-8;
systemPrms.amp=13e-12;
systemPrms.ks=20.0;
systemPrms.v=(4/3)*pi*0.7^3*1e-15;

u0 = systemPrms.rest;
tspan = [0,0.01];
dt = 0.000001;
tvec = tspan(1):dt:tspan(2);
NT = length(tvec);

% Instant Solution
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',dt);
instantSols = zeros(NT, NA);
for na = 1:NA
    systemPrms.amp = amp(na);
    isol = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms),tspan,u0,odeOptions);
    instantSols(:,na) = deval(isol,tvec);
end

plot(tvec,instantSols);


%% For range of [B]t's, optimize kon and koff for 3 ica amplitudes

NBT = 20;
bt = logspace(-5,-1,NBT);
konRange = [5e7 5e12];
kdRange = [5e-6 3e-3];

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






%% Functions Needed for Optimization

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
function dy = massActionDynamics(t,y,systemPrms,prms) %#ok
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

% Gives solution with mass-action buffer but Fixed BT
function dy = massActionDynamics_FixBT(t,y,systemPrms,bt,prms)
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
    association = y(1)*(bt-y(2))*prms(1);
    dissociation = y(2)*(prms(1)*prms(2)); % Prms 2 is Kd
    bufferExchange = association - dissociation;
    
    % Derivative
    dy(1,1) = curr - extrusion - bufferExchange;
    dy(2,1) = bufferExchange;
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








