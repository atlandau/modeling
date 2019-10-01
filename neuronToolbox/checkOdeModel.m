
[odeModel,cellMorph,cellStructure,compStructure] = generateCell();
NC = size(odeModel.physPrm,1);
iState = odeModel.physPrm(:,2);

tolerance = 1e-6;
odeOptions = odeset('AbsTol',tolerance,'RelTol',tolerance);

tspan = [0 0.2];
odeProblem = @(t,y) odeCellFunction(t,y,odeModel);
[t,v] = ode23s(odeProblem,tspan,iState,odeOptions);



%
idxSpines = find(cellStructure.type==2);
idxSpineCompartment = cell2mat(compStructure.gci(idxSpines));
idxParent = cell2mat(compStructure.parent(idxSpines));

figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.3 0.64 0.60 0.34]);

subplot(1,2,1); 
hold on;
plot(1000*t,1000*v(:,idxSpineCompartment(1)),'color','r','linewidth',1.5);
plot(1000*t,1000*v(:,idxParent(1)),'color','k','linewidth',1.5);
title(compStructure.name{idxSpines(1)});

subplot(1,2,2); 
hold on;
plot(1000*t,1000*v(:,idxSpineCompartment(2)),'color','r','linewidth',1.5);
plot(1000*t,1000*v(:,idxParent(2)),'color','k','linewidth',1.5);
title(compStructure.name{idxSpines(2)});


%%
imagesc((0:NC-1)*dx,1000*t,1000*v)
colormap('hot');
colorbar
