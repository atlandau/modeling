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




dv = (p.vc(t)/p.rs - vm/p.rs - vm/p.rm + p.em/p.rm)/p.cm;




