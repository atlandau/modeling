

hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/portFromJulia';
fpath = '/Users/landauland/Documents/Julia/calciumBuffering/latexDiscussion/images';

% Spine Parameters
systemPrms.beta=1400.0;
systemPrms.rest=5.0e-8;
systemPrms.amp=6e-12;
systemPrms.ks=20.0;
systemPrms.v=(4/3)*pi*0.7^3*1e-15;

u0 = systemPrms.rest;
tspan = [0.001,0.01];
dt = 0.000001;
tvec = tspan(1):dt:tspan(2);

% Instant Solution
odeOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'MaxStep',dt);
[trueTime,trueSol] = ode23s(@(t,y) instantBufferDynamics(t,y,systemPrms),tspan,u0,odeOptions);
plot(trueTime,trueSol);


%% Grid Search
NKON = 10;
NKD = 25;
NBT = 20;
kon = logspace(5,10,NKON);
kd = logspace(-8,-2,NKD);
bt = logspace(-6,-1,NBT);

costs = zeros(NKON,NKD,NBT);
trajectories = zeros(length(trueTime),NKON,NKD,NBT);
bufferTraj = zeros(length(trueTime),NKON,NKD,NBT);

fprintf('working... ');
msg = '';
for nkon=1:NKON
    for nkd=1:NKD
        fprintf(repmat('\b',1,length(msg)));
        msg = sprintf('Kon: %d/%d, Kd: %d/%d...\n',nkon,NKON,nkd,NKD);
        fprintf(msg);
        
        for nbt=1:NBT
            koff = kon(nkon)*kd(nkd);
            cprms = [kon(nkon) koff bt(nbt)];
            kappa = bt(nbt) / (systemPrms.rest + kd(nkd));
            cu0 = [systemPrms.rest; kappa*systemPrms.rest]; % Don't let bad initialization penalize parameters
            cSolution = ode23s(@(t,y)massActionDynamics(t,y,systemPrms,cprms), trueTime, cu0, odeOptions);
            ctraj = deval(cSolution,trueTime);
            trajectories(:,nkon,nkd,nbt) = ctraj(1,:);
            bufferTraj(:,nkon,nkd,nbt) = ctraj(2,:);
            costs(nkon,nkd,nbt) = sum((ctraj(1,:)' - trueSol).^2);
        end
    end
end
fprintf('finished.\n\n');

% save('gridSearch_190813.mat','trajectories','bufferTraj','costs','bk','ks','-v7.3');



%% plot 5 best
NBest = 100;
[~,idxAscending] = sort(costs(:));
idxBest = idxAscending(1:NBest);
[konBest,kdBest,btBest] = ind2sub(size(costs),idxBest);

[konAll,kdAll,btAll] = ind2sub(size(costs),idxAscending);
kdOrdered = kd(kdAll);
btOrdered = bt(btAll);

fontSize = 16;
rbMap = cat(2, linspace(0,1,NBest)',zeros(NBest,1),linspace(1,0,NBest)');

figure(3); clf; 
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.7 0.5]);

subplot(2,3,1); 
hold on;
plot(trueTime,trueSol,'linewidth',1.5,'color','k');
for nbest = 1:NBest
    plot(1e3*trueTime,1e6*trajectories(:,konBest(nbest),kdBest(nbest),btBest(nbest)),'linewidth',1.5,'color',rbMap(nbest,:));
end
xlabel('Time (ms)');
ylabel('[Ca] µM');
title(sprintf('Instant & %d Best Traj''s',NBest));
set(gca,'fontsize',fontSize);

subplot(2,3,4); 
hold on;
for nbest = 1:NBest
    cerror = trajectories(:,konBest(nbest),kdBest(nbest),btBest(nbest)) - trueSol(:);
    plot(1e3*trueTime,1e6*cerror,'linewidth',1.5,'color',rbMap(nbest,:));
end
ylim([-0.1 0.3]*1e-1);
xlabel('Time (ms)');
ylabel('Error [Ca] µM');
title(sprintf('Error %d Best Traj''s',NBest));
set(gca,'fontsize',fontSize);

NSPlots = numel(costs)/20;
xLim = [1 numel(costs)];
yLim = [1e-15 1e-6];
subplot(2,3,[2 5]); hold on;
p = patch([1 1 NSPlots NSPlots],yLim([1 2 2 1]),'y');
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';
plot(costs(idxAscending),'color','k','linewidth',1.5); 
xlim(xLim); 
ylim(yLim); 
set(gca,'xscale','log'); 
set(gca,'yscale','log'); 
xlabel('Solution Number'); 
ylabel('Cost (sum(MSE))');
title('Cost');
set(gca,'fontsize',fontSize);

textToPlot = {
    'patch indicates number of solutions';
    '   shown in scatter plot on right'};
text(0.5,0.1,textToPlot,'units','normalized');

subplot(2,3,[3 6]); 
mkrSize = 50*norman(flipud(costs(idxAscending(1:NSPlots))))+1;
mkrColor = hot(NSPlots);
% pp = plot(kdOrdered,btOrdered,'linestyle','none','marker','o','color','k','markerfacecolor','k','markersize',3);
ss = scatter(kdOrdered(1:NSPlots),btOrdered(1:NSPlots),mkrSize,mkrColor,'filled');
alpha(gca,0.1);
set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('Kd');
ylabel('[B]t');
title('Kd vs. Bt');
set(gca,'fontsize',fontSize);


inset = axes('Position',[0.43 0.6 0.1 0.25]);
box on;
% hold on;
% p = patch([1 1 NSPlots NSPlots],yLim([1 2 2 1]),'y');
% p.FaceAlpha = 0.1;
% p.EdgeColor = 'none';
plot(costs(idxAscending),'color','k','linewidth',1.5);
xlim(xLim);
ylim(yLim);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
% set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('Cost w/ linear scaling');
set(gca,'fontsize',fontSize);

% print(gcf,'-painters',fullfile(hpath,'costReport_bigGridSearch'),'-djpeg');


%% Functions Needed for Optimization

function out = ica(t,amp)
    % function giving calcium current with fixed temporal properties and
    % user-defined amplitude
    current = amp*exp(-(t-2e-3).^2 / (0.55e-3)^2);
    offset = amp*exp(-(1e-3-2e-3).^2 / (0.55e-3)^2);
    adjustedCurrent = (current - offset) .* (t>=1e-3).*(t<=3e-3);
    out = adjustedCurrent .* (adjustedCurrent>=0);
end

% Gives solution with instant buffer
function dy = instantBufferDynamics(t,y,prm)
    % Current/Extrusion Term
    restCurrent = prm.beta/prm.ks * prm.rest;
    curr = restCurrent + ica(t,prm.amp)/(2*96485*prm.v*prm.ks);
    extrusion = prm.beta/prm.ks * y(1);
    % Derivative
    dy = curr - extrusion;
end

% Gives solution with mass-action buffer
function dy = massActionDynamics(t,y,systemPrms,prms)
    % Parameters: [v,ks,beta,rest,amp]
    % Kon,Koff,Bt are optimization variables
    
    % y(1) = cafree
    % y(2) = boundBuffer
    % y(3) = secondBoundBuffer
    
    % Current/Extrusion Term
    restCurrent = systemPrms.beta * systemPrms.rest;
    curr = restCurrent + ica(t,systemPrms.amp)/(2*96485*systemPrms.v);
    extrusion = systemPrms.beta * y(1);
    
    % Association/Dissociation with quasi-instant buffer
    association = y(1)*(prms(3)-y(2))*prms(1);
    dissociation = y(2)*prms(2);
    bufferExchange = association - dissociation;
    
    % Derivative
    dy(1,1) = curr - extrusion - bufferExchange;
    dy(2,1) = bufferExchange;
end

% Gives solution with mass-action buffer - if no saturation is allowed
function dy = massActionDynamics_noSat(t,y,systemPrms,prms)
    % Parameters: [v,ks,beta,rest,amp]
    % Kon,Koff,Bt are optimization variables
    
    % y(1) = cafree
    % y(2) = boundBuffer
    % y(3) = secondBoundBuffer
    
    % Current/Extrusion Term
    restCurrent = systemPrms.beta * systemPrms.rest;
    curr = restCurrent + ica(t,systemPrms.amp)/(2*96485*systemPrms.v);
    extrusion = systemPrms.beta * y(1);
    
    % Association/Dissociation with quasi-instant buffer
    association = y(1)*prms(3)*prms(1);
    dissociation = y(2)*prms(2);
    bufferExchange = association - dissociation;
    
    % Derivative
    dy(1,1) = curr - extrusion - bufferExchange;
    dy(2,1) = bufferExchange;
end


function cost = objFcn(prms,systemPrms,odeOptions,trueTime,trueSol,tSpan,u0)
    odeSol = ode23s(@(t,y) massActionDynamics(t,y,systemPrms,prms), tSpan,u0,odeOptions);
    testSol = deval(odeSol, trueTime);
    cost = testSol(1,:)' - trueSol;
end








