
prm.rneck = 1000e6;
prm.rhead = 2e12 / (4*pi*0.5^2); % Ohm-µm^2 / surfaceArea
prm.chead = 1e-14 * (4*pi*0.5^2);
prm.rdend = sqrt(2e12*1e6/2)/(2*pi*(5/2)^(3/2));
prm.cdend = 3e-12; %kind of a guess
prm.e = -0.07;

gpeak = 1e-9;
tpeak = 0.002;
prm.gsyn = @(t) (t>0).* gpeak .*exp(1)/tpeak.*t.*exp(-t./tpeak);

gSyn = @(t,tPeak,gPeak) (t>0) .* gPeak .* exp(1)/tPeak .* t .* exp(-t./tPeak); % Conductance function (t in milliseconds)
tspan = [0 0.03];
tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);
iState = [prm.e; prm.e];

odeProblem = @(t,v) vcSpineDendrite(t,v,prm);
[t,v] = ode45(odeProblem,tspan,iState,odeOptions);
plot(1000*t,1000*v)
