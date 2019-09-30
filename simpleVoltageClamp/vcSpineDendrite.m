function dv = vcSpineDendrite(t,v,prm)
% dv = vcSpineDendrite(t,v,prm)
%
% t is time in s
% v is the membrane potential in volts
% prm is the parameter structure
%   - prm.rs = spine resistance
%   - prm.cs = spine capacitance
%   - prm.es = spine reversal
%   - prm.rn = spineneck resistance
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
