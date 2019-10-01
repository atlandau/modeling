
clear cellMorphology

cellMorphology.soma.morphid = 0;
cellMorphology.soma.type = 0;
cellMorphology.soma.diameter = 25;

cellMorphology.amain.morphid = 1;
cellMorphology.amain.type = 1;
cellMorphology.amain.diameter = 5;
cellMorphology.amain.length = 20000;
cellMorphology.amain.link = 0;
cellMorphology.amain.location = 0;

cellMorphology.cc.morphid = 2;
cellMorphology.cc.type = 12;
cellMorphology.cc.link = 1;
cellMorphology.cc.location = 10000;
cellMorphology.cc.ccCommand = @(t) 100e-12;


dx = 50;
odeCellModel = generateCell(cellMorphology,[],dx);
NC = size(odeCellModel.physPrm,1);
iState = odeCellModel.physPrm(:,2);

tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);

tspan = [0 0.1];
odeProblem = @(t,y) odeCellFunction(t,y,odeCellModel);
[t,v] = ode45(odeProblem,tspan,iState,odeOptions);


figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.3 0.64 0.60 0.34]);

subplot(1,3,1);
imagesc((0:NC-1)*dx,1000*t,1000*v);
colormap('hot');
colorbar;

subplot(1,3,2);
plotvc((0:NC-1)*dx,1000*v(1:100:end,:),gcf);
xlabel('Distance (µm)');
ylabel('Potential');

subplot(1,3,3);
plotvc(1000*t,1000*v(:,1:10:end),gcf);




%% Single Spine Synapse

clear cellMorphology

cellMorphology.soma.morphid = 0;
cellMorphology.soma.type = 0;
cellMorphology.soma.diameter = 25;

cellMorphology.amain.morphid = 1;
cellMorphology.amain.type = 1;
cellMorphology.amain.diameter = 5;
cellMorphology.amain.length = 20000;
cellMorphology.amain.link = 0;
cellMorphology.amain.location = 0;

cellMorphology.spine1.morphid = 3;
cellMorphology.spine1.type = 2;
cellMorphology.spine1.diameter = 0.8;
cellMorphology.spine1.rneck = 500e6;
cellMorphology.spine1.link = 1;
cellMorphology.spine1.location = 10000;

cellMorphology.synapse1.morphid = 4;
cellMorphology.synapse1.type = 3;
cellMorphology.synapse1.link = 3;
cellMorphology.synapse1.properties = [0, 2e-9, 0.001, 0];

profile on;

dx = 50;
odeCellModel = generateCell(cellMorphology,[],dx);
NC = size(odeCellModel.physPrm,1);
iState = odeCellModel.physPrm(:,2);

tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);

tspan = [0 0.01];
odeProblem = @(t,y) odeCellFunction(t,y,odeCellModel);
[t,v] = ode45(odeProblem,tspan,iState,odeOptions);

profile viewer;


figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.3 0.64 0.60 0.34]);

subplot(1,3,1); hold on;
plot(1000*t,1000*v(:,402),'color','r','linewidth',1.5);
plot(1000*t,1000*v(:,201),'color','k','linewidth',1.5);
legend('Spine','Dendrite','location','northeast');

subplot(1,3,2);
plotvc((0:NC-2)*dx,1000*v(1:100:end,1:end-1),gcf);
xlim([0 (NC-2)*dx]);
xlabel('Distance (µm)');
ylabel('Potential');

subplot(1,3,3);
plotvc(1000*t,1000*v(:,1:10:end),gcf);





























