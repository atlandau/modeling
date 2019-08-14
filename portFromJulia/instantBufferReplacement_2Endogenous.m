



% Spine Parameters
systemPrms.beta=1400.0;
systemPrms.rest=5.0e-8;
systemPrms.amp=30e-12;
systemPrms.ks=20.0;
systemPrms.v=(4/3)*pi*0.7^3*1e-15;

u0 = systemPrms.rest;
tspan = [0.001,0.005];
dt = 0.000001;
tvec = tspan(1):dt:tspan(2);

odeOptions = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'MaxStep',dt);
[trueTime,trueSol] = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms),tspan,u0,odeOptions);

z0 = [systemPrms.rest; 10*systemPrms.rest; 10*systemPrms.rest];
objective = @(prms) objFcn(prms,systemPrms,odeOptions,trueTime,trueSol,tspan,z0);
iPrm = [5e9;  875000; 3e-3; 5e7; 125000; 50e-6];
lowerBound = [5e6; 1000; 1e-6; 5e6; 1000; 1e-6];
upperBound = [5e20; 9e15; 1; 5e20; 9e15; 1e-3];

lsqOptions = optimoptions('lsqnonlin','OptimalityTolerance',1e-20,'FunctionTolerance',1e-20,'StepTolerance',1e-20);

sol = lsqnonlin(objective, iPrm, lowerBound, upperBound, lsqOptions);

bestSolution = ode23s(@(t,y)massActionDynamics(t,y,systemPrms,sol), tSpan, z0);
bestSol = deval(bestSolution(1,:), trueTime);


figure(1); 
clf;
hold on;

plot(trueTime,trueSol,'linewidth',1.5,'color','k');
plot(trueTime, bestSol(1,:), 'linewidth',1.5,'color','r');

fprintf('Buffer1: Kon=%.2fx10^9, Kd=%.1fµM, [B]t=%.3fmM\n',sol(1)*1e-9,sol(2)/sol(1)*1e6, sol(3)*1e3);
fprintf('Buffer2: Kon=%.2fx10^9, Kd=%.1fµM, [B]t=%.3fmM\n',sol(4)*1e-9,sol(5)/sol(4)*1e6, sol(6)*1e3);



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
    association1 = y(1)*(prms(3)-y(2))*prms(1);
    dissociation1 = y(2)*prms(2);
    bufferExchange1 = association1 - dissociation1;
    
    association2 = y(1)*(prms(6)-y(3))*prms(4);
    dissociation2 = y(3)*prms(5); 
    bufferExchange2 = association2 - dissociation2;
    
    % Derivative
    dy(1,1) = curr - extrusion - bufferExchange1 - bufferExchange2;
    dy(2,1) = bufferExchange1;
    dy(3,1) = bufferExchange2;
end

function cost = objFcn(prms,systemPrms,odeOptions,trueTime,trueSol,tSpan,u0)
    odeSol = ode23s(@(t,y) massActionDynamics(t,y,systemPrms,prms), tSpan,u0,odeOptions);
    testSol = deval(odeSol, trueTime);
    cost = testSol(1,:)' - trueSol;
end









