pkFluor = 0.6; 
gmax = 3;
kdFluo5f = 800e-9;
totalFluor = 300e-6; % M
caExtracellular = 1.5; % mM
caProportionNMDAR = 1 / (1 + 65.6/(0.5*11.3*caExtracellular)); % From Jahr/Stevens 1993
F = 96485; % C/mole

spineArea = 4*pi*0.75^2; % (in µM^2);
spineVolume = 4/3*pi*0.75^3; % (in µM^3)
specificCapacitance = 1e-6*(1/10000)^2; % (in F/µM^2)
capSpine = specificCapacitance * spineArea; % in F

rfluor = pkFluor / gmax;
caConc = rfluor*totalFluor;
naConc = caConc * (1/caProportionNMDAR);
totalChargeDensity = caConc*2 + naConc; % Mole / L
totalCharge = totalChargeDensity * F * spineVolume*1e-15;

voltage = totalCharge / capSpine

