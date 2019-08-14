function dv = vcSomaDendriteDEQ(t,v,prm)
% dv = vcdiffeq(t,vm,p)
%
% t is time in ms
% vm is the membrane potential in volts
% prm is the parameter structure
%   - prm.rp = access
%   - prm.rs = soma resistance
%   - prm.cs = soma capacitance
%   - prm.rd = dendrite resistance
%   - prm.cd = dendrite capacitance
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
%
% -- the equations --
% Irp = Irm + Icm + Ira        |    KCL
% Ira = Ird + Icd
% 
% dVs/dt = (Vh-Vs)/CsRp - Vs/CsRs - (Vs-Vd)/CsRa
% dVd/dt = (Vs-Vd)/CdRa - Vd/CdRd

% assume that Es and Ed are 0...

if t>=0.005
%     disp()
end

dvs = ((prm.vc(t) - v(1))/prm.rp - v(1)/prm.rs - (v(1)-v(2)/prm.ra))/prm.cs; % soma 
dvd = ((v(1)-v(2))/prm.ra - v(2)/prm.rd)/prm.cd; % dendrite

dv = [dvs; dvd]; 

if any(isnan(dv))
%     disp()
end



