function dydt = bufferDynamics_011(t,y,ica,p)
% dydt = bufferDynamics_011(t,y,ica,p)
% 
% Nomenclature (011): [endogenousBuffer,fluorescentBuffer,singleSpine]
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
%   p.beta is the rate (1/ms) of the yoked extrusion mechanisms
%   p.rest is the [Ca]rest - required to offset input current


% Derivatives
currentTerm = p.beta*p.rest + ica(t,p.amplitude)/(96485*p.v); % change due to current
associationTerm = p.kon * y(1) * (p.bconc - y(2)); % change due to association with calcium and buffer
dissociationTerm = p.koff * y(2); % change due to dissociation of calcium-bound buffer
extrusionTerm = p.beta * y(1); % change due to extrusion pumps

% Final Output
dydt = [currentTerm - associationTerm + dissociationTerm - extrusionTerm; associationTerm - dissociationTerm];



%NOTE: this is somewhat optimized. vectorization could probably be better
