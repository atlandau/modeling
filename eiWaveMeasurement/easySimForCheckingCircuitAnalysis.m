
% Simulation Parameters
dt = 0.001e-3; % s
ds = 1;
T = 50e-3; % s
tspan = [0 T];

% System Parameters
p.rs = 10e6; % access resistance (Ohms) 
p.rm = 500e6; % input resistance (Ohms)
p.cm = 20e-12; % cell capacitance (F)
p.em = -70e-3; % rest potential

% Set up voltage clamp hold potential
vcType = 'step'; % options: {'sine','step'}
vcModDepth = 0e-3; % volts
vcModPeriod = 2e-3; % seconds
vcModCenter = -35e-3; % center voltage
p.vc = @(t) vcModCenter + vcModDepth/2 * sin(2*pi*t/vcModPeriod);

% Run Simulation
iState = p.vc(0); % start at 0mV (this model assumes 0mV as rest potential) 
[t,vm] = eulerapp(@(t,vm) vcdiffeq(t,vm,p),tspan,iState,dt,ds); % do euler approximation 

Ivc = (p.vc(0:dt*ds:T)' - vm)/p.rs; % Voltage Clamp (Measured) Current 
Ir = vm/p.rm; % Membrane Current 
Ic = p.cm * diff(vm)/dt; % Capacitive Current 


% Compute Estimate from parameters
if strcmp(vcType,'sine')
    % Homogeneous Solution
    tauTotal = (p.rm*p.rs)/(p.rm+p.rs) * p.cm;
    
    emrpTerm = p.em*(p.rs/(p.rs+p.rm));
    responseHomogeneous = emrpTerm - (emrpTerm - iState)*exp(-t/tauTotal);
    
    % Step Solution
    emv0term = (p.em*p.rs + vcModCenter*p.rm) / (p.rm+p.rs);
    responseStep = emv0term - (emv0term-iState)*exp(-t/tauTotal);
    
    % Sine Solution    
    w = 2*pi/vcModPeriod;
    estVm = @(t,center,depth,period) center + depth/2 * sin(2*pi*t/period);
    
    realDamping = (p.rm + p.rs)/(p.rm);
    imagDamping = w * p.rs * p.cm;    
    magnitudeResponse = vcModDepth / sqrt(realDamping^2 + imagDamping^2);
    timeDelay = atan(imagDamping/realDamping)/(2*pi) * vcModPeriod;    
    
    responseSine = estVm(t-timeDelay,0,magnitudeResponse,vcModPeriod);
    
    estResponse = responseStep + responseSine;% + responseHomogeneous;
end

%
g = figure(1);
clf;
set(gcf,'units','normalized','outerposition',[0.3 0.1 0.4 0.7]);

subplot(2,1,1);
hold on;
plot(1e3*t,1e3*p.vc(t),'color','k','linewidth',1.5);
plot(1e3*t,1e3*vm,'color','r','linewidth',1.5);
plot(1e3*t,1e3*estResponse,'color','c','linewidth',1.5);
xlabel('Time (ms)');
ylabel('V_m (mV)');
title('Voltage');
legend('V_{hold}','V_{m}','V_{est}','location','northeast');
set(gca,'fontsize',16);

subplot(2,1,2);
hold on;
plot(1e3*t,1e3*responseHomogeneous,'k','linewidth',1.5);
plot(1e3*t,1e3*responseStep,'r','linewidth',1.5);
plot(1e3*t,1e3*responseSine,'c','linewidth',1.5);
xlabel('Time(ms)');
ylabel('V_m');
title('Estimate of Each Solution');
legend('Homogeneous','Step','Sine','linewidth',1.5);
set(gca,'fontsize',16);

% 
% [r,lags] = xcorr(norman(Ivc),norman(vm));
% bestLag = lags(r==max(r));
% sampleDelay = -bestLag;
% timeDelay = t(abs(sampleDelay));
% angleDelay = -sign(bestLag)*timeDelay/vcModPeriod*360;
% fprintf('angle delay: %d degrees\n',round(angleDelay));


%% - plot result -
f = figure(1);
clf;
set(gcf,'units','normalized','outerposition',[0.3 0.1 0.4 0.7]);

subplot(3,1,1);
hold on;  
plot(1e3*t,1e3*p.vc(t),'color','k','linewidth',1.5);
plot(1e3*t,1e3*vm,'color','r','linewidth',1.5);
if strcmp(vcType,'sine')
    plot(1e3*t,1e3*estimateResponse,'color','c','linewidth',1.5);
end
xlabel('Time (ms)');
ylabel('V_m (mV)');
title('Voltage');
if strcmp(vcType,'sine')
    legend('V_{hold}','V_{m}','V_{est}','location','northeast');
else
    legend('V_{hold}','V_{m}','location','northeast');
end
set(gca,'fontsize',16);

subplot(3,1,2);
hold on;
plot(1e3*t,1e9*Ivc,'color','k','linewidth',1.5);
plot(1e3*t,1e9*Ir,'color','r','linewidth',1.5);
plot(1e3*t(1:end-1),1e9*Ic,'color','c','linewidth',1.5);
xlabel('Time (ms)');
ylabel('Current (nA)');
title('Currents');
legend('I_{vc}','I_{rm}','I_{cm}','location','northeast');
set(gca,'fontsize',16);

subplot(3,1,3);
hold on;
plot(1e3*t,norman(Ivc),'color','k','linewidth',1.5);
plot(1e3*t,norman(vm),'color','r','linewidth',1.5);
xlabel('Time (ms)');
ylabel('Normalized');
title('Phase Delay I_{vc}, V_m');
legend('I_{vc}','V_m','location','northeast');
set(gca,'fontsize',16);










































