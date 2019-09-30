function defaultMembrane = membraneParameterDefaults(type)
% defaultMembrane = membraneParameterDefaults(type)
% 
% returns structure containing default parameters for membrane physiology
% type determines if it's a passive or active membrane, default is passive 

if nargin==0, type = 0; end

switch type
    case {0,'passive'}        
        defaultMembrane.sMemRes = 2e12; % Ohm-µm^2
        defaultMembrane.eMem = -70e-3; % V
        defaultMembrane.sMemCap = 1e-14; % F/µm^2
        defaultMembrane.axialRes = 1e6; % Ohm-µm
end
