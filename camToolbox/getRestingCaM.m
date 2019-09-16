function [camDistribution,sol] = getRestingCaM(camConcentration,restCalcium)

camProperties = getCaMProps();

tspan = [0 10];
if nargin==3, tspan = [0 endSimulation]; end
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);
iState = [restCalcium;camConcentration;zeros(8,1)];

odeProblem = @(t,y) camODE_baselineEstimation(t,y,camProperties);
sol = ode23s(odeProblem,tspan,iState,odeOptions);
camDistribution = sol.y(2:end,end);

