function dydt = bufferDynamics_KK1(t,y,ica,p)
% dydt = bufferDynamics_KK1(t,y,ica,p)
% 
% Nomenclature (KK1): [InstantaneousEndogenousBuffer,InstantaneousFluorescentBuffer,singleSpine]
% 
% t is time (which is irrelevant because this process is memoryless)
% y is the data 
%   y(1): [Ca]free -- concentration free calcium in cell
%   y(2): [CaB] -- concentration of calcium bound to buffer
%
% ica an inline function giving calcium current at time t
%
% p provides parameters required for computation
%   p.v is volume of compartment
%   p.kon is association rate of buffer
%   p.koff is dissociation rate of buffer
%   p.bconc is the total concentration of the buffer
%   p.gamma is the rate (1/ms) of the yoked extrusion mechanisms
%   p.rest is the [Ca]rest - required to offset input current
%   p.ks is the instantaneous endogenous buffer capacity


kappa = p.bconc / (p.koff/p.kon + p.rest); % assume constant kappa derived from resting value

% Derivatives
newCurrent = ica(t,p.amplitude,p.ks)/(96485*p.v);

currentTerm = p.gamma/(1+p.ks)*p.rest + ; % change due to current
associationTerm = p.kon * y(1) * (p.bconc - y(2)); % change due to association with calcium and buffer
dissociationTerm = p.koff * y(2); % change due to dissociation of calcium-bound buffer
extrusionTerm = p.gamma/(1+p.ks) * y(1); % change due to extrusion pumps

% Final Output
dydt = [currentTerm - associationTerm + dissociationTerm - extrusionTerm; associationTerm - dissociationTerm];



%NOTE: this is somewhat optimized. vectorization could probably be better


