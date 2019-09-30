

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

tspan = [0 0.05];

vrest = [sdprm.es sdprm.ed];
odeOptions = odeset('AbsTol',1e-6,'RelTol',1e-6);
odeProblem = @(t,v) vcSomaDendrite(t,v,sdprm);
[t,v] = ode23s(odeProblem,tspan,vrest,odeOptions);
curr = (sdprm.vc(t) - v(:,1)) / sdprm.rp;

figure(1); clf;
subplot(3,1,1);
plot(1e3*t,1e3*v,'linewidth',1.5);
ylabel('Membrane Potential (mV)');
set(gca,'fontsize',16);

subplot(3,1,2);
plot(1e3*t,1e12*curr,'color','r','linewidth',1.5);
xlabel('Time (ms)');
ylabel('VC Current (pA)');
set(gca,'fontsize',16);

subplot(3,1,3);
plot(1e3*t,1e3*sdprm.vc(t),'color','r','linewidth',1.5);
xlabel('Time (ms)');
ylabel('V Command');
set(gca,'fontsize',16);

%% Soma + Dendrite
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
    cdvs = (prm.vc(t) - v(1))/prm.rp - (v(1) - prm.es)/prm.rs - (v(1)-v(2))/prm.ra;
    cdvd = (v(1) - v(2))/prm.ra - (v(2) - prm.ed)/prm.rd;
    dv = [cdvs/prm.cs; cdvd/prm.cs]; 
end






