function dy = instantBufferDynamics(t,y,prm)
% Parameters: [v,ks,beta,rest,amp,kon,koff]
% y = [Ca]

% Current/Extrusion Term
restCurrent = prm.beta/(prm.ks+1) * prm.rest;
curr = restCurrent + ica(t,prm.amp)/(2*96485*prm.v*(prm.ks+1));
extrusion = prm.beta/(prm.ks+1) * y(1);

% Derivative
dy = curr - extrusion;
