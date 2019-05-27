function dydt = bufferDynamics_121(t,y,ica,p)
% dydt = bufferDynamics_121(t,y,ica,p)
% 
% Nomenclature (121): [endogenousBuffer,fluorescentBuffer,singleSpine]
%                                       2: Simple two state GCaMP model
% 
% t is time (which is irrelevant because this process is memoryless)
% y is the data 
%   y(1): [Ca]free -- concentration free calcium in cell
%   y(2): [CaB] -- concentration of calcium bound to buffer
%   y(3): [CaG] -- concentration of calcium bound to GCaMP --
%                  note that this is cooperative with 4 binding sites
%
% ica an inline function giving calcium current at time t
%
% p provides parameters required for computation
%   p.v is volume of compartment
%   p.kon is association rate of buffer
%   p.koff is dissociation rate of buffer
%   p.bconc is the total concentration of the buffer
%   p.beta is the rate (1/ms) of the yoked extrusion mechanisms


% -------------------------------------------------------------------------
% this is somewhat optimized. the vectorization could probably be better


% if ica(t)>0
%     disp('here')
% end

cafree = y(1);
endogenousBound = y(2);
gcampBound = y(3);

% Derivatives
currentTerm = p.beta*5e-8 + ica(t)/(96485*p.v); % change due to current
extrusionTerm = p.beta * cafree; % extrusion

% Endogenous Reactions
endogenousTerm = p.endogenousKOn * cafree * (p.endConc - endogenousBound) - p.endogenousKOff * endogenousBound;

% GCaMP Reactions -- 4x for calcium, 1x for GCaMP
gcampTerm = p.gcampKOn * cafree^p.gcampHill * (p.gcampConc - gcampBound) - p.gcampKOff * gcampBound; 


% Final Output
dydt = [
    currentTerm - extrusionTerm - endogenousTerm - 4*gcampTerm;
    extrusionTerm;
    gcampTerm];


