
%before running IV curve, check the input resistance
clear cellMorphology
cellMorphology = cellMorphologyDefaults();

cellMorphology.cc.morphid = 13;
cellMorphology.cc.type = 12;
cellMorphology.cc.link = 0;
cellMorphology.cc.ccCommand = @(t) 100e-12;

clear defaultMembrane
defaultMembrane.sMemRes = 0.25e12; % Ohm-µm^2
defaultMembrane.eMem = -70e-3; % V
defaultMembrane.sMemCap = 1e-14; % F/µm^2
defaultMembrane.axialRes = 0.7e6; % Ohm-µm

tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);

tspan = [0 0.05];
dx = 50;
odeCellModel = generateCell(cellMorphology,defaultMembrane,dx);

iState = odeCellModel.physPrm(:,2);
    
odeProblem = @(t,y) odeCellFunction(t,y,odeCellModel);
[t,v] = ode45(odeProblem,tspan,iState,odeOptions);
plot(1000*t,1000*v(:,1),'k');
fprintf(1,'Input Resistance: %.2f\n',1e-6*(v(end,1)-v(1,1))/100e-12);

    
    
%% Here- we vary the voltage clamp holding potential and record the different currents

clear cellMorphology
cellMorphology = cellMorphologyDefaults();

cellMorphology.synapse1.morphid = 13;
cellMorphology.synapse1.type = 3;
cellMorphology.synapse1.link = 11;
cellMorphology.synapse1.properties = [0, 1e-9, 0.001, 0.05];

cellMorphology.vc.morphid = 14;
cellMorphology.vc.type = 11;
cellMorphology.vc.link = 0;
cellMorphology.vc.vcAccess = 10e6; 
cellMorphology.vc.vcCommand = @(t) -0.07;

clear defaultMembrane
defaultMembrane.sMemRes = 0.25e12; % Ohm-µm^2
defaultMembrane.eMem = -70e-3; % V
defaultMembrane.sMemCap = 1e-14; % F/µm^2
defaultMembrane.axialRes = 0.7e6; % Ohm-µm


tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);

tspan = [0 0.1];
dx = 50;
[odeCellModel,cellMorph,cellStructure,compStructure] = generateCell(cellMorphology,[],dx);
NC = size(odeCellModel.physPrm,1);
iState = odeCellModel.physPrm(:,2);

somaIdx = 1;
dendIdx = 9;
spineIdx = 12;


holdingPotentials = -0.08:0.02:-0.04;
NH = length(holdingPotentials);
t = cell(1,NH);
v = cell(1,NH);
vcCurr = cell(1,NH);
syCurr = cell(1,NH);
for nh = 1:NH
    fprintf(1,'Holding Potential %d/%d...\n',nh,NH);
    cellMorphology.vc.vcCommand = @(t) holdingPotentials(nh);
    [odeCellModel,cellMorph,cellStructure,compStructure] = generateCell(cellMorphology,[],dx);
    
    odeProblem = @(t,y) odeCellFunction(t,y,odeCellModel);
    [t{nh},v{nh}] = ode45(odeProblem,tspan,iState,odeOptions);
    
    vcCurr{nh} = zeros(numel(t{nh}),1);
    syCurr{nh} = zeros(numel(t{nh}),1);
    for i = 1:numel(t{nh})
        vcCurr{nh}(i) = (odeCellModel.voltageClamp{1}(t{nh}(i))-v{nh}(i,somaIdx))/cellMorphology.vc.vcAccess;
        syCurr{nh}(i) = odeCellModel.synConductance{1}(t{nh}(i)).*(odeCellModel.synReversal{1}-v{nh}(i,spineIdx)); 
    end
end
fprintf(1,'done!!!\n');

%%
somaMap = repmat(linspace(0.75,0,NH)',1,3);
spineMap = [linspace(0.75,0,NH)',zeros(NH,2)];
somaMap = spineMap;

figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.3 0.39 0.42 0.59]);

subplot(2,2,1); 
hold on;
for nh = 1:NH
    plot(1000*t{nh},1000*v{nh}(:,somaIdx),'color',somaMap(nh,:),'linewidth',1.5);
end
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Soma Voltage');
set(gca,'fontsize',16);

subplot(2,2,3); 
hold on;
for nh = 1:NH
    plot(1000*t{nh},1000*v{nh}(:,spineIdx),'color',spineMap(nh,:),'linewidth',1.5);
end
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Spine Voltage');
set(gca,'fontsize',16);

subplot(2,2,2); 
hold on;
for nh = 1:NH
    tBaseline = find(t{nh}>=0.049,1);
    plot(1000*t{nh},1e12*(vcCurr{nh}-vcCurr{nh}(tBaseline)),'color',somaMap(nh,:),'linewidth',1.5);
end
xlim([40 80]);
xlabel('Time (ms)');
ylabel('Current (pA)');
title('VoltageClamp Measurement');
set(gca,'fontsize',16);

subplot(2,2,4); 
hold on;
for nh = 1:NH
    tBaseline = find(t{nh}>=0.049,1);
    plot(1000*t{nh},1e12*(syCurr{nh}-syCurr{nh}(tBaseline)),'color',spineMap(nh,:),'linewidth',1.5);
end
xlim([40 80]);
xlabel('Time (ms)');
ylabel('Current (pA)');
title('Synaptic Current Measurement');
set(gca,'fontsize',16);
















































