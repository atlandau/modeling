
cellMorphology.soma.morphid = 0;
cellMorphology.soma.type = 0;
cellMorphology.soma.diameter = 25;
cellMorphology.amain.morphid = 1;
cellMorphology.amain.type = 1;
cellMorphology.amain.diameter = 5;
cellMorphology.amain.length = 10000;
cellMorphology.amain.link = 0;
cellMorphology.amain.location = 0;

dx = 50;
odeCellModel = generateCell(cellMorphology,[],dx);
NC = size(odeCellModel.physPrm,1);
iState = odeCellModel.physPrm(:,2);

tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);

tspan = [0 0.2];
odeProblem = @(t,y) odeCellFunction(t,y,odeCellModel);
[t,v] = ode23s(odeProblem,tspan,iState,odeOptions);


figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.3 0.64 0.60 0.34]);

subplot(1,3,1);
imagesc(1:NC,1000*t,1000*v);
colormap('hot');
colorbar;

subplot(1,3,2);
plotvc((0:NC-1)*10,1000*v(1:10:end,:),gcf);
xlabel('Distance (µm)');
ylabel('Potential');

subplot(1,3,3);
plotvc(1000*t,1000*v(:,1:10:end),gcf);



%%
imagesc((0:NC-1)*dx,1000*t,1000*v)
colormap('hot');
colorbar
