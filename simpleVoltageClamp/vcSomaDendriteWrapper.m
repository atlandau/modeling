
% I set prm.ra to inf so the dendrite is effectively not connected. If you
% want to connect it set it to something like 50 MOhm


% to calculate access, divide the initial current amplitude after the step
% by the amplitude of the voltage change. As you rp/rs -> 0, it will get
% more and more accurate.

% Pipette resistance
prm.rp = 10e6;

% Soma resistance, capacitance, reversal
prm.rs = 100e6;
prm.cs = 150e-12; 
prm.es = -70e-3; 

% Axial resistance between soma and dendrite
prm.ra = 100e6; %50e6; 

% Dendrite resistance, capacitance, reversal
prm.rd = 100e6; 
prm.cd = 100e-12;
prm.ed = -70e-3;

% Voltage-Clamp Parameters
holdBaseline = -70e-3; 
holdStep = -30e-3;
startVStep = 50e-3; % start voltage step 
endVStep = 100e-3; % end voltage step 
rcAmplitude = holdStep - holdBaseline;

% inline function describing hold potential as a function of time
prm.vc = @(t) holdBaseline + (t>=startVStep & t<endVStep)*rcAmplitude;

% Setup diffeq simulation 
tspan = [0 0.2]; % seconds
iState = [prm.es; prm.ed];
tolerance = 1e-8; % how precise to do the simulation
odeOptions = odeset('AbsTol',1e-6,'RelTol',1e-6);

% ode45 requires the function to take two inputs, just time and value
% this parameterizes the system with 'prm' so we can feed this into ode45
odeProblem = @(t,v) vcSomaDendrite(t,v,prm); 
[t,v] = ode45(odeProblem,tspan,iState,odeOptions);

% Voltage Clamp Current
curr = (prm.vc(t) - v(:,1)) / prm.rp;

% Estimate Access Resistance
% (note- I'm doing this with a short window because it's more realistic-
% digital samples don't take 0 time but in the simulation the first step
% will effectively happen instantaneously)
tStep = find(t>=startVStep,1,'first'); % first sample during voltage step
tAfter = find(t>=startVStep+0.00025,1,'first'); % some µs after voltage step starts
accessEstimate = (holdStep-holdBaseline)/mean(curr(tStep:tAfter)); % estimate of access resistance in MOhm
perError = 100*(accessEstimate - prm.rp)/prm.rp;

% Estimate Input Resistance
lastSample = find(t<endVStep,1,'last'); % Can use a single point because there's no noise in the simulation
ssCurrent = curr(lastSample);
inputResistance = rcAmplitude/(ssCurrent-curr(1)) - accessEstimate;

% Capacitance estimate
ft = fittype('a*exp(-x/tau)');
fo = fitoptions(ft);
fo.StartPoint = [1 0.01]; 
fo.Upper = [2 1];
fo.Lower = [0 0];
decayTime = t(tStep:lastSample) - t(tStep); % Start at 0
decayValue = (curr(tStep:lastSample) - ssCurrent) / (curr(tStep)-ssCurrent);
fObj = fit(decayTime(:),decayValue(:),ft,fo);
tau = fObj.tau;
conventionalCapEstimate = tau / accessEstimate; % What people usually do
moreAccurateCapEstimate = tau * (1/accessEstimate + 1/inputResistance); % Pipette and soma resistance combined

trueInputResistance = 1/( (1/prm.rs) + (1/(prm.ra+prm.rd)));
fprintf(1,'Estimate of access resistance: %.1f M%c (True Value: %.1f M%c)\n',1e-6*accessEstimate,937,1e-6*prm.rp,937);
fprintf(1,'Estimate of input resistance: %.1f M%c (True Value: %.1f M%c)\n',1e-6*inputResistance,937,1e-6*trueInputResistance,937);
fprintf(1,'Estimate of capacitance: %.1f pF (True Value: %.1f pF)\n',1e12*moreAccurateCapEstimate,1e12*prm.cs);
fprintf(1,'Bad Estimate of capacitance: %.1f pF (True Value: %.1f pF)\n',1e12*conventionalCapEstimate,1e12*prm.cs);
fprintf(1,'\n');


% Plotting
figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.68 0.28 0.31 0.55]);

subplot(3,1,1); hold on;
plot(1e3*t,1e3*v(:,1),'color','k','linewidth',1.5);
plot(1e3*t,1e3*v(:,2),'color','r','linewidth',1.5);
ylabel('mV');
title('Membrane Potential');
legend('Soma','Dendrite','location','northeast');
set(gca,'fontsize',16);

subplot(3,1,2); hold on;
plot(1e3*t,1e12*curr,'color','k','linewidth',1.5);
xlabel('Time (ms)');
ylabel('pA');
title('VoltageClamp Current');
set(gca,'fontsize',16);

subplot(3,1,3);
plot(1e3*t,1e3*prm.vc(t),'color','k','linewidth',1.5);
ylim([-75 5]);
xlabel('Time (ms)');
ylabel('mV');
title('Voltage Command');
set(gca,'fontsize',16);



%% ------ 
function dv = vcSomaDendrite(t,v,prm)
    % t is time in ms
    % vm is the membrane potential in volts
    % prm is the parameter structure
    %   - prm.rp = access
    %   - prm.rs = soma resistance
    %   - prm.cs = soma capacitance
    %   - prm.es = soma reversal
    %   - prm.rd = dendrite resistance
    %   - prm.cd = dendrite capacitance
    %   - prm.ed = dendrite reversal
    %   - prm.vc = inline for time-dependent change
    %
    % Differential equations describing voltage clamp circuit
    %            -
    %           ---  
    %            |
    %           Vvc
    %            | == Vh
    %            Rp
    %            |
    % Vs == -----------------------Ra----------------------- == Vd
    %       |            |                    |            |
    %       Rs           |                    Rd           |
    %       |            Cs                   |            Cd
    %       Es           |                    Ed           |
    %       |            |                    |            |
    %       --------------                    -------------- 
    %            |                                   |
    %           ---                                 ---
    %            -                                   -
    % Irp = Irs + Ics + Ira
    % Ira = Ird + Icd
    % Cs * dVs/dt = (Vh-Vs)/Rp - (Vs-Es)/Rs - (Vs-Vd)/Ra
    % Cd * dVd/dt = (Vs-Vd)/Ra - (Vd-Ed)/Rd
    holdVoltage = prm.vc(t);
    somaVoltage = v(1);
    dendVoltage = v(2);
    cdvs = (holdVoltage-somaVoltage)/prm.rp - (somaVoltage-prm.es)/prm.rs - (somaVoltage-dendVoltage)/prm.ra; % soma 
    cdvd = (somaVoltage-dendVoltage)/prm.ra - (dendVoltage-prm.ed)/prm.rd; % dendrite
    dv = [cdvs/prm.cs; cdvd/prm.cd]; 
end

















