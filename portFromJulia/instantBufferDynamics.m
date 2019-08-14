function dy = instantBufferDynamics(t,y,prm)
% Parameters: [v,ks,beta,rest,amp,kon,koff]
% y = [Ca]

% Current/Extrusion Term
restCurrent = prm.beta/prm.ks * prm.rest;
curr = restCurrent + ica(t,prm.amp)/(96485*prm.v*prm.ks);
extrusion = prm.beta/prm.ks * y(1);

% Derivative
dy = curr - extrusion;
