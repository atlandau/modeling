function dv = vcdiffeq(t,vm,p)
% dv = vcdiffeq(t,vm,p)
%
% super stupid model- assumes Erest is 0V
%
% t is time in ms
% vm is the membrane potential in volts
% p is the parameter structure
%   - p.rs = access
%   - p.rm = input resistance
%   - p.cm = capacitance
%   - p.vc = scalar for constant value or inline for time-dependent change


% Set up voltage clamp potential
if isa(p.vc, 'function_handle')
    vc = p.vc(t);
else
    vc = p.vc;
end

% So fucking simple HA
dv = (vc/p.rs - vm/p.rs - vm/p.rm)/p.cm;



