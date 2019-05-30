

p.rs = 10e6; % access resistance (Ohms) 
p.rm = 150e6; % input resistance (Ohms)
p.cm = 100e-12; % cell capacitance (F)

% Voltage Clamp step at 5ms
vcStepTime = 5e-3; % s
vcHold = 100e-3; % Volts
p.vc = @(t) vcHold * (t>=vcStepTime); % heaviside


dt = 0.01e-3; % s
ds = 1;
T = 20e-3; % s
tspan = [0 T];
iState = 0;

[t,vm] = eulerapp(@(t,vm) vcdiffeq(t,vm,p),tspan,iState,dt,ds); % do euler approximation

Ivc = (p.vc(0:dt*ds:T)' - vm)/p.rs; % Voltage Clamp Current
Ir = vm/p.rm; % Membrane Current
Ic = p.cm * diff(vm)/dt; % Capacitive Current


% - plot result -
f = figure;
clf;

subplot(2,1,1);
plot(1e3*t,1e3*vm,'color','k','linewidth',1.5)
xlabel('Time (ms)');
ylabel('V_m (mV)');
title('Cell Voltage');
set(gca,'fontsize',16);

subplot(2,1,2);
hold on;
plot(1e3*t,1e9*Ivc,'color','k','linewidth',1.5);
plot(1e3*t,1e9*Ir,'color','r','linewidth',1.5);
plot(1e3*t(1:end-1),1e9*Ic,'color','c','linewidth',1.5);
xlabel('Time (ms)');
ylabel('Current (nA)');
title('Currents');
legend('I_{vc}','I_{rm}','I_{cm}','location','northeast');
set(gca,'fontsize',16);



%% -- functions that the above code needs to run --
function [t,y,dy] = eulerapp(ode,tspan,initState,dt,ds)
    % [t,y,dy] = eulerapp(ode,tspan,initState,dt,ds)
    % 
    % t - time vector
    % y - values
    % dy - derivatives
    %
    % ode is an inline function 
    % tspan is the start and end point (will truncate if not multiple of dt)
    % initState is initial state of function
    % dt is time step used
    % ds is a downsample factor for having a highres time step but not stupidly
    %    large data sizes - works simply (1:ds:end)

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

%{
% Differential equation describing voltage clamp circuit
% ______________________
%            -
%           ---  
%            |
%           Ivc
%            |        ---------------- Vvc
%            Rs
%            |
%    --------------------------------- Vm
%    |               |
%    Rm              Cm
%    |               |
%    -----------------
%            |
%           ---
%            -
% ______________________
%
%
% -- the equations --
% Irs = Irm + Icm
% Irs = (Vvc - Vm) / Rs
% Irm = Vm/Rm
% Icm = Cm * dVm/dt
%
% dVm/dt = (Vvc/Rs) - Vm/Rs - Vm/Rm
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
    %   - p.vc = scalar for constant value or inline for time-dependent change


    % Set up voltage clamp potential
    if isa(p.vc, 'function_handle')
        vc = p.vc(t);
    else
        vc = p.vc;
    end

    % ODE
    dv = (vc/p.rs - vm/p.rs - vm/p.rm)/p.cm;
end






