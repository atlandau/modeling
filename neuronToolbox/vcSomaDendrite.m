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
%   - prm.vc = inline equation for time-dependent change of holding voltage
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
% -- the equations --
% Irp = Irm + Icm + Ira        |    KCL
% Ira = Ird + Icd
% dVs/dt = (Vh-Vs)/CsRp - Vs/CsRs - (Vs-Vd)/CsRa
% dVd/dt = (Vs-Vd)/CdRa - Vd/CdRd

cdvs = (prm.vc(t) - v(1))/prm.rp - (v(1) - prm.es)/prm.rs - (v(1)-v(2))/prm.ra;
cdvd = (v(1) - v(2))/prm.ra - (v(2) - prm.ed)/prm.rd;
dv = [cdvs/prm.cs; cdvd/prm.cs]; 
