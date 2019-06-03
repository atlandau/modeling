
fpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/eiWaveMeasurement/figures_analysis_190602';

% Simulation Parameters
dt = 0.0001e-3; % s
ds = 1;
T = 100e-3; % s
tspan = [0 T];

% System Parameters
p.rs = 10e6; % access resistance (Ohms) 
p.rm = 150e6; % input resistance (Ohms)
p.cm = 100e-12; % cell capacitance (F)
p.em = -70e-3; % rest potential

% Set up voltage clamp hold potential
vcType = 'sine'; % options: {'sine','step'}
vcModDepth = 30e-3; % volts
vcModPeriod = 2e-3; % seconds
vcModCenter = -70e-3; % center voltage
p.vc = @(t) vcModCenter + vcModDepth/2 * sin(2*pi*t/vcModPeriod);

plotCycles = 4;
vcModPeriod = 2e-3;%[0.5 1 2 5 10]*1e-3;
NP = length(vcModPeriod);
response=cell(NP,3);
for np = 1:NP

    p.vc = @(t) vcModCenter + vcModDepth/2 * sin(2*pi*t/vcModPeriod(np));
    
    % Run Simulation
    iState = p.vc(0); % start at 0mV (this model assumes 0mV as rest potential) 
    [t,vm] = eulerapp(@(t,vm) vcdiffeq(t,vm,p),tspan,iState,dt,ds); % do euler approximation 

    Ivc = (p.vc(0:dt*ds:T)' - vm)/p.rs; % Voltage Clamp (Measured) Current 
    Ir = vm/p.rm; % Membrane Current 
    Ic = p.cm * diff(vm)/dt; % Capacitive Current 


    % Compute Estimate from parameters
    if strcmp(vcType,'sine')
        tauTotal = (p.rm*p.rs)/(p.rm+p.rs) * p.cm;

        % Step Solution
        emv0term = (p.em*p.rs + vcModCenter*p.rm) / (p.rm+p.rs);
        responseStep = emv0term - (emv0term-iState)*exp(-t/tauTotal);

        % Sine Solution    
        w = 2*pi/vcModPeriod(np);
        estVm = @(t,center,depth,period) center + depth/2 * sin(2*pi*t/period);

        realDamping = (p.rm + p.rs)/(p.rm);
        imagDamping = w * p.rs * p.cm;    
        magnitudeResponse = vcModDepth / sqrt(realDamping^2 + imagDamping^2);
        timeDelay = atan(imagDamping/realDamping)/(2*pi) * vcModPeriod(np);    

        responseSine = estVm(t-timeDelay,0,magnitudeResponse,vcModPeriod(np));

        estResponse = responseStep + responseSine;% + responseHomogeneous;
    end
    
    nSamples = round(plotCycles*vcModPeriod(np)/dt);
    response{np,1} = p.vc(t(1:nSamples))';
    response{np,2} = vm(1:nSamples);
    response{np,3} = estResponse(1:nSamples)';
end

%%
plotIdx = vcModPeriod==2e-3;
plotNumSamples = length(response{plotIdx,1});
plotTime = dt:dt:dt*plotNumSamples;

ff = figure(1);
clf;
set(gcf,'units','normalized','outerposition',[0.3 0.1 0.3 0.7]);

subplot(2,1,1);
hold on;
plot(1e3*plotTime,1e3*response{plotIdx,1},'color','k','linewidth',1.5);
plot(1e3*plotTime,1e3*response{plotIdx,2},'color','r','linewidth',1.5);
plot(1e3*plotTime,1e3*response{plotIdx,3},'color','c','linewidth',1.5);
xlabel('Time (ms)');
ylabel('V_m (mV)');
title('Voltage Response and Estimate to 500Hz Input');
legend('V_{hold}','V_{m}','V_{est}','location','northeast');
set(gca,'fontsize',16);

cmap = varycolor(NP);
subplot(2,1,2);
hold on;
for np = 1:NP
    cNumSamples = length(response{np,1});
    ds = plotCycles/cNumSamples;
    tt = ds:ds:plotCycles;
    plot(tt,1e3*(response{np,3}-response{np,2}),'color',cmap(np,:),'linewidth',1.5);
end
set(gca,'xtick',0:1:plotCycles);
xlabel('Cycles');
ylabel('Error (mV)');
title('Error in Estimates');
legend(cellfun(@(c) sprintf('%4s ms/cyc',num2str(round(c*1e4)/10)),num2cell(vcModPeriod),'uni',0),'location','southeast');
set(gca,'fontsize',16);

% print(gcf,'-painters',fullfile(fpath,'ResponseToSineInput'),'-djpeg');


%%
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
    plot(1e3*t,1e3*estResponse,'color','c','linewidth',1.5);
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










































