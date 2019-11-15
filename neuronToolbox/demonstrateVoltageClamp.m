
% measure effect of voltage clamp of soma on compartmentalized EPSC


clear cellMorphology

cellMorphology.soma.morphid = 0;
cellMorphology.soma.type = 0;
cellMorphology.soma.diameter = 12;

cellMorphology.synapse1.morphid = 1;
cellMorphology.synapse1.type = 3;
cellMorphology.synapse1.link = 0;
cellMorphology.synapse1.properties = [0, 1e-9, 0.001, 0];

cellMorphology.vc.morphid = 2;
cellMorphology.vc.type = 11;
cellMorphology.vc.link = 0;
cellMorphology.vc.vcAccess = 10e6; 
cellMorphology.vc.vcCommand = @(t) -0.070;

dx = 50;
[odeCellModel,cellMorph,cellStructure,compStructure] = generateCell(cellMorphology,[],dx);
NC = size(odeCellModel.physPrm,1);
iState = odeCellModel.physPrm(:,2);

tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);

tspan = [0 0.02];
odeProblem = @(t,y) odeCellFunction(t,y,odeCellModel);
[t,v] = ode45(odeProblem,tspan,iState,odeOptions);
vcCurr = (odeCellModel.voltageClamp{1}(t)-v(:,1))/cellMorphology.vc.vcAccess;
syCurr = zeros(numel(t),1);
for i = 1:numel(t), syCurr(i) = odeCellModel.synConductance{1}(t(i)).*(odeCellModel.synReversal{1}-v(i,1)); end

%
figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.3 0.64 0.60 0.34]);

subplot(2,1,1);
hold on;
plot(1000*t,1000*v,'color','k','linewidth',1.5);

subplot(2,1,2);
hold on;
plot(1000*t,1e12*syCurr,'color','b','linewidth',1.5);
plot(1000*t,-1e12*vcCurr,'color','k','linewidth',1.5);


% 
% subplot(1,3,2);
% plotvc((0:NC-2)*dx,1000*v(1:100:end,1:end-1),gcf);
% xlim([0 (NC-2)*dx]);
% xlabel('Distance (µm)');
% ylabel('Potential');
% 
% subplot(1,3,3);
% plotvc(1000*t,1000*v(:,1:10:end),gcf);
% 




























