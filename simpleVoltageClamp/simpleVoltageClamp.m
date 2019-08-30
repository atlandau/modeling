
% Describes the voltage command
baselineHold = -70e-3;
rcAmplitude = 10e-3; % V
rcStart = 5e-3; % s
rcEnd = 25e-3; % s
voltageCommand = @(t) baselineHold + rcAmplitude.*(t>=rcStart).*(t<=rcEnd); % inline function giving voltage step

% Parameters
sprm.rp = 10e6; % pipette resistance
sprm.rm = 100e6; % membrane resistance
sprm.em = -75e-3; % GHK
sprm.cm = 150e-12; % membrane capacitance
sprm.vc = voltageCommand;

tspan = [0 0.1]; 

vrest = (baselineHold*sprm.rm + sprm.em*sprm.rp)/(sprm.rp + sprm.rm); % steady-state offset
odeOptions = odeset('AbsTol',1e-6,'RelTol',1e-6);
odeProblem = @(t,v) vcSoma(t,v,sprm);
[t,v] = ode45(odeProblem,tspan,vrest,odeOptions);
dv = vcSoma(t,v,sprm);
curr = (sprm.vc(t) - v) / sprm.rp; % pipette current (see Irp below...)

figure(1); clf;
subplot(2,1,1);
plot(1e3*t,1e3*v,'k');
ylabel('Membrane Potential (mV)');
set(gca,'fontsize',16);

subplot(2,1,2);
plot(1e3*t,1e12*curr,'r');
xlabel('Time (ms)');
ylabel('VC Current (pA)');
set(gca,'fontsize',16);




%% ODEs
function dv = vcSoma(t,vm,p)
    % dv = vcSoma(t,vm,p)
    %
    % t is time in ms
    % vm is the membrane potential in volts
    % p is the parameter structure
    %   - p.rp = access resistance
    %   - p.rm = input resistance
    %   - p.cm = capacitance
    %   - p.em = rest potential
    %   - p.vc = inline for time-dependent change
    %
    % Differential equation describing voltage clamp circuit
    % 
    %            -
    %           ---  
    %            |
    %           Ivc
    %            |  ------ - ---------------- Vc
    %            Rp
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
    % -- the equations --
    % Irp = Irm + Icm         |    KCL
    %
    % Irp = (Vc - Vm) / Rp    |    Ivc = Irs
    % Irm = (Vm-Em)/Rm
    % Icm = Cm * dVm/dt
    %  
    % dVm/dt = (Vc/Rp - Vm/Rp - Vm/Rm + Em/Rm)/Cm
    dv = (p.vc(t)/p.rp - vm/p.rp - vm/p.rm + p.em/p.rm)/p.cm;
end





