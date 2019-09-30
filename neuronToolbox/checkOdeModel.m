

[odeCellModel,cellMorph,cellTable,cStructure] = generateCell();
NC = size(odeCellModel.physPrm,1);
iState = odeCellModel.physPrm(:,2);

tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);

tspan = [0 0.2];
odeProblem = @(t,y) odeCellFunction(t,y,odeCellModel);
[t,v] = ode23s(odeProblem,tspan,iState,odeOptions);


%%
idxSpines = cellTable(2,:)==2;
idxSpineCompartment = cell2mat(cStructure(3,idxSpines));
idxParent = cell2mat(cStructure(4,idxSpines));

figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.3 0.64 0.60 0.34]);

%subplot(1,3,1); 
hold on;
plot(1000*t,1000*v(:,idxSpineCompartment),'color','k','linewidth',1.5);
plot(1000*t,1000*v(:,idxParent),'color','r','linewidth',1.5);



%%
imagesc((0:NC-1)*dx,1000*t,1000*v)
colormap('hot');
colorbar
