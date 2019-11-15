function dv = vcSomaDendrite(t,v,prm)
% t is time in ms
% vm is the membrane potential in volts
% prm is the parameter structure
%   - prm.rp = access
%   - prm.rs = soma resistance
%   - prm.cs = soma capacitance
%   - prm.es = soma reversal
%   - prm.rd = dendrite resistance
%   - prm.cd = dendrite capacitance
%   - prm.ed = dendrite reversal
%   - prm.vc = inline for time-dependent change
%
% Differential equations describing voltage clamp circuit
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
% Irp = Irs + Ics + Ira
% Ira = Ird + Icd
% Cs * dVs/dt = (Vh-Vs)/Rp - (Vs-Es)/Rs - (Vs-Vd)/Ra
% Cd * dVd/dt = (Vs-Vd)/Ra - (Vd-Ed)/Rd
holdVoltage = prm.vc(t);
somaVoltage = v(1);
dendVoltage = v(2);
cdvs = (holdVoltage-somaVoltage)/prm.rp - (somaVoltage-prm.es)/prm.rs - (somaVoltage-dendVoltage)/prm.ra; % soma 
cdvd = (somaVoltage-dendVoltage)/prm.ra - (dendVoltage-prm.ed)/prm.rd; % dendrite
dv = [cdvs/prm.cs; cdvd/prm.cd]; 