
sprm.rp = 10e6;
sprm.rm = 100e6;
sprm.em = -70e-3;
sprm.cm = 150e-12;
sprm.vc = @(t) sprm.em + 10e-3.*(t>=5e-3).*(t<=25e-3); % inline function giving voltage step

tspan = [0 0.1];

vrest = sprm.em;
odeOptions = odeset('AbsTol',1e-8,'RelTol',1e-8);
odeProblem = @(t,v) vcSoma(t,v,sprm);
[t,v] = ode23s(odeProblem,tspan,vrest,odeOptions);
dv = vcSoma(t,v,sprm);
curr = (sprm.vc(t) - v) / sprm.rp;

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


%% Soma + Dendrite
sdprm.rp = 10e6;
sdprm.rs = 80e6;
sdprm.cs = 100e-12;
sdprm.es = -70e-3;
sdprm.ra = 50e6;
sdprm.rd = 80e6;
sdprm.cd = 100e-12;
sdprm.ed = -70e-3;
sdprm.vc = @(t) sdprm.es + 10e-3.*(t>=5e-3).*(t<=25e-3); % inline function giving voltage step

vrest = [sdprm.es sdprm.ed];
odeOptions = odeset('AbsTol',1e-10,'RelTol',1e-10);
odeProblem = @(t,v) vcSomaDendrite(t,v,sdprm);
[t,v] = ode23s(odeProblem,tspan,vrest,odeOptions);
curr = (sdprm.vc(t) - v(1)) / sdprm.rp;

figure(1); clf;
subplot(2,1,1);
plot(1e3*t,1e3*v);
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
    % Irs = Irm + Icm         |    KCL
    %
    % Irs = (Vc - Vm) / Rp    |    Ivc = Irs
    % Irm = (Vm-Em)/Rm
    % Icm = Cm * dVm/dt
    %  
    % dVm/dt = (Vc/Rp - Vm/Rs - Vm/Rm + Em/Rm)/Cm
    dv = (p.vc(t)/p.rp - vm/p.rp - vm/p.rm + p.em/p.rm)/p.cm;
end

function dv = vcSomaDendrite(t,v,prm)
    % dv = vcSomaDendrite(t,v,prm)
    %
    % t is time in s
    % v is the membrane potential in volts
    % prm is the parameter structure
    %   - prm.rp = access resistance
    %   - prm.rs = soma resistance
    %   - prm.cs = soma capacitance
    %   - prm.es = soma reversal
    %   - prm.ra = axial resistance
    %   - prm.rd = dendrite resistance
    %   - prm.cd = dendrite capacitance
    %   - prm.ed = dendrite reversal
    %   - prm.vc = inline for time-dependent change
    %
    % Differential equations describing voltage clamp circuit
    % 
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
    %
    % -- the equations --
    % Irp = Irm + Icm + Ira        |    KCL
    % Ira = Ird + Icd
    % 
    % dVs/dt = (Vh-Vs)/CsRp - Vs/CsRs - (Vs-Vd)/CsRa
    % dVd/dt = (Vs-Vd)/CdRa - Vd/CdRd
    dvs = ((v(1) - prm.vc(t))/prm.rp + (prm.es - v(1))/prm.rs + (v(2)-v(1))/prm.ra)/prm.cs;
    dvd = ((v(1) - v(2))/prm.ra + (prm.ed - v(2))/prm.rd)/prm.cd;
    dv = [dvs; dvd]; 
end






