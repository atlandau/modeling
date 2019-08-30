
prm.rp = 10e6;
prm.rs = 300e6;
prm.cs = 30e-12;
prm.rd = 100e6;
prm.ra = 200e6;
prm.cd = 90e-12;

vstart = 0;
vhold = 10e-3;
vtime = 0;
prm.vc = @(t) (t>=vtime)*vhold; % step to vhold at vtime

dt = 0.000001e-3; % seconds (very short...)
ds = 1000;
T = 0.001; % seconds
tspan = [0 T];

iState = [0;0];
[t,v,dv] = eulerapp(@(t,v) vcSomaDendriteDEQ(t,v,prm),tspan,iState,dt,ds,1);

ms = 1e3*t;

figure(1); clf;
hold on;
plot(ms,1e3*v(:,1),'k','linewidth',1.5);
plot(ms,1e3*v(:,2),'r','linewidth',1.5);
xlabel('Time (ms)');
ylabel('mV');
title('Response to Voltage Step');
legend('Soma','Dendrite','location','northeast');
set(gca,'fontsize',16);


odeProblem = @(t,v) vcdiffeq(t,v,prm);
[t,v] = ode23s(odeProblem,[0 1],0);






%% ODEs

function dv = vcdiffeq(t,vm,p)
    % dv = vcdiffeq(t,vm,p)
    %
    % t is time in ms
    % vm is the membrane potential in volts
    % p is the parameter structure
    %   - p.rs = access
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
    % -- the equations --
    % Irs = Irm + Icm         |    KCL
    %
    % Irs = (Vc - Vm) / Rs    |    Ivc = Irs
    % Irm = (Vm-Em)/Rm
    % Icm = Cm * dVm/dt
    %  
    % dVm/dt = (Vc/Rs - Vm/Rs - Vm/Rm + Em/Rm)/Cm
    dv = (p.vc(t)/p.rs - vm/p.rs - vm/p.rm + p.em/p.rm)/p.cm;
end





