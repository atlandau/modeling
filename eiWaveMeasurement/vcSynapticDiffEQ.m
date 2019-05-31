function dv = vcSynapticDiffEQ(t,vm,p,syn)
% dv = vcSynapticDiffEQ(t,vm,p,syn)
%
% t is time in ms
% vm is the membrane potential in volts
% p is the parameter structure for membrane properties
%   - p.rs = access
%   - p.rm = input resistance
%   - p.cm = capacitance
%   - p.em = rest potential
%   - p.vc = inline for time-dependent change
%
% syn is the parameter structure for synaptic properties
%   - syn.ee = excitatory reversal
%   - syn.ei = inhibitory reversal
%   - syn.tvec = time vector of corresponding conductances in syn.ge/gi
%   - syn.ge = excitatory conductance as function of tvec
%   - syn.gi = inhibitory conductance as function of tvec
%   - *** -
%     synaptic conductances are sometimes stochastic so there isn't an
%     obvious way to pass an inline function through, this method uses a
%     time vector and a conductance array 
%
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
%    |    |     |    |   
%    |    Rm    Ge   Gi
%    Cm   |     |    |
%    |    Em    Ee   Ei
%    |    |     |    |
%    ----------------- - ---------------- Ground
%            |
%           ---
%            -
%
%
% -- the equations --
% Irs = Irm + Icm + Ige + Igi         |    KCL
%
% Irs = (Vc - Vm) / Rs                |    Ivc = Irs
% Irm = (Vm-Em)/Rm
% Icm = Cm * dVm/dt
% Ige = (Vm-Ee)*Ge
% Igi = (Vm-Ei)*Gi
% 
% Cm * dVm/dt = (Vc-Vm)/Rs - (Vm-Em)/Rm - (Vm-Ee)*Ge - (Vm-Ei)*Gi;

if numel(syn.tvec) ~= numel(syn.ge) || numel(syn.tvec) ~= numel(syn.gi)
    error('Synaptic time vector (syn.tvec) has to be same size as synaptic conductances (syn.ge, syn.gi)');
end

tIndex = find(syn.tvec >= t,1);
cge = syn.ge(tIndex);
cgi = syn.gi(tIndex);

Irs = (p.vc(t) - vm)/p.rs;
Irm = (vm-p.em)/p.rm;
Ige = (vm-syn.ee)*cge;
Igi = (vm-syn.ei)*cgi;

dv =  (Irs - Irm - Ige - Igi) / p.cm;



