

% System Parameters
p.rs = 10e6; % access resistance (Ohms) 
p.rm = 300e6; % input resistance (Ohms)
p.cm = 200e-12; % cell capacitance (F)
p.em = -70e-3; % rest potential

% Set up voltage clamp hold potential
vcType = 'sine'; % options: {'sine','step'}
switch vcType
    case 'sine'
        vcModDepth = 10e-3; % volts
        vcModPeriod = 0.5e-3; % seconds
        vcModCenter = -70e-3; % center voltage
        p.vc = @(t) vcModCenter + vcModDepth/2 * sin(2*pi*t/vcModPeriod);
    case 'step'
        vcStepTime = 5e-3; % s
        vcHold = -70e-3; % Volts
        p.vc = @(t) vcHold * (t>=vcStepTime); % heaviside
end

% Simulation Parameters
dt = 0.01e-3; % s
ds = 1;
T = 20e-3; % s
tspan = [0 T];
iState = p.vc(0);%-0.070774287780911; % start at 0mV (this model assumes 0mV as rest potential) 

[t,vm] = eulerapp(@(t,vm) vcdiffeq(t,vm,p),tspan,iState,dt,ds); % do euler approximation 

% Compute Currents - (see circuit model below)
Ivc = (p.vc(0:dt*ds:T)' - vm)/p.rs; % Voltage Clamp (Measured) Current 
Ir = vm/p.rm; % Membrane Current 
Ic = p.cm * diff(vm)/dt; % Capacitive Current 

% Compute Estimate from parameters
if strcmp(vcType,'sine')
    w = 2*pi/vcModPeriod;
    estVm = @(t,center,depth,period) center + depth/2 * sin(2*pi*t/period);
    
    realDamping = (p.rm + p.rs)/(p.rm);
    imagDamping = w * p.rs * p.cm;
    
    dampFactor = sqrt(realDamping^2 + imagDamping^2);
    magnitudeResponse = vcModDepth / dampFactor;
    angleResponse = atan(imagDamping/realDamping); 
    timeDelay = angleResponse/(2*pi)*vcModPeriod;
    estimateResponse = estVm(t-timeDelay,vcModCenter,magnitudeResponse,vcModPeriod);
end
[r,lags] = xcorr(norman(Ivc),norman(vm));
bestLag = lags(r==max(r));
sampleDelay = -bestLag;
timeDelay = t(abs(sampleDelay));
angleDelay = -sign(bestLag)*timeDelay/vcModPeriod*360;
fprintf('angle delay: %d degrees\n',round(angleDelay));


% - plot result -
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




%% -- functions that the above code needs to run --

%{
% Differential equation describing voltage clamp circuit
% 
%            -
%           ---  
%            |
%           Ivc
%            |  ------ - ---------------- Vc
%            Rs
%            |
%    ----------------- - ---------------- Vm
%    |               |
%    Rm              |
%    |               Cm
%    Em              |
%    |               |
%    ----------------- - ---------------- Ground
%            |
%           ---
%            -
%
%
% -- the equations --
% Irs = Irm + Icm         |    KCL
%
% Irs = (Vc - Vm) / Rs    |    Ivc = Irs
% Irm = (Vm-Em)/Rm
% Icm = Cm * dVm/dt
% 
% 
% dVm/dt = (Vc/Rs - Vm/Rs - Vm/Rm + Em/Rm)/Cm
%
%}
function dv = vcdiffeq(t,vm,p)
    % dv = vcdiffeq(t,vm,p)
    %
    % super stupid model- assumes Erest is 0V
    %
    % t is time in ms
    % vm is the membrane potential in volts
    % p is the parameter structure
    %   - p.rs = access
    %   - p.rm = input resistance
    %   - p.cm = capacitance
    %   - p.em = rest potential
    %   - p.vc = inline for time-dependent change
    dv = (p.vc(t)/p.rs - vm/p.rs - vm/p.rm + p.em/p.rm)/p.cm;
end


% numerical approximation wrapper
function [t,y,dy] = eulerapp(ode,tspan,initState,dt,ds)
    % [t,y,dy] = eulerapp(ode,tspan,initState,dt,ds)
    % 
    % t - time vector
    % y - values
    % dy - derivatives
    %
    % ode is an inline function 
    % tspan is the start/end point
    % initState is initial state of function
    % dt - time step used
    % ds - downsample factor for having a highres time step but not
    %      stupidly large data sizes - works simply (1:ds:end)

    if nargin<5, ds = 1; end
    if rem(tspan(2),dt*ds)~=0, error('time vector must be chosen perfectly'); end

    t = tspan(1):dt*ds:tspan(2); % time vector
    NT = length(t); % number of data points
    NV = length(initState); % number of values in function
    
    % preallocate
    y = zeros(NT,NV);
    dy = zeros(NT-1,NV);

    % setup initial state
    initState = initState(:)';
    y(1,:) = initState;
    tempy = initState;

    % Report progress to screen
    fprintf(1,'| euclidean approximation working... ');
    msg = '';
    counter = 0;
    percentage = 10 - 9*(NT*ds>=1e7); % 10% if fast, 1% if slow

    % loop through and perform euclidean approximation
    for i = 1:NT-1, ctime = t(i);

        % Progress report
        if 100*i/NT > counter
            fprintf(1,repmat('\b',1,length(msg)-1));
            msg = sprintf('%d%%%%',counter);
            fprintf(1,msg);
            counter = counter+percentage;
        end

        % Loop through each sub-time point
        for j = 1:ds, computeTime = ctime + dt*(j-1);
            tempdy = ode(computeTime,tempy);
            tempy = tempy + tempdy(:)'*dt;
        end
        y(i+1,:) = tempy;
        dy(i,:) = tempdy;
    end

    fprintf(1,repmat('\b',1,length(msg)-1));
    fprintf(1,' finished.\n');
end



