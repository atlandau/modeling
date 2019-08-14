

hpath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/portFromJulia';
fpath = '/Users/landauland/Documents/Julia/calciumBuffering/latexDiscussion/images';

% Spine Parameters
systemPrms.beta=1400.0;
systemPrms.rest=5.0e-8;
systemPrms.amp=3.3e-12;
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
NBK = 50;
NKS = 21;
bk = logspace(5,16,NBK);
ks = linspace(15,25,NKS);

costs = zeros(NBK,NKS);
trajectories = zeros(length(trueTime),NBK,NKS);
bufferTraj = zeros(length(trueTime),NBK,NKS);

fprintf('working... ');
msg = '';
for nbk= 1:NBK
    for nks=1:NKS
        fprintf(repmat('\b',1,length(msg)));
        msg = sprintf('BK: %d/%d, KS: %d/%d...',nbk,NBK,nks,NKS);
        fprintf(msg);
        
        cprms = [bk(nbk) ks(nks)];
        cu0 = [systemPrms.rest; ks(nks)*systemPrms.rest]; % Don't let bad initialization penalize parameters
        cSolution = ode23s(@(t,y)massActionDynamicsRPRM(t,y,systemPrms,cprms), trueTime, cu0, odeOptions);
        ctraj = deval(cSolution,trueTime);
        trajectories(:,nbk,nks) = ctraj(1,:);
        bufferTraj(:,nbk,nks) = ctraj(2,:);
        costs(nbk,nks) = sum((ctraj(1,:)' - trueSol).^2);
    end
end
fprintf('finished.\n\n');

% save('gridSearchRPRM_190809.mat','trajectories','bufferTraj','costs','bk','ks','-v7.3');


%% Visualize Results

rbMap = @(N) cat(2, linspace(1,0,N)',zeros(N,1), linspace(0,1,N)');
fontSize = 26;

figure(1); clf; 
set(gcf,'units','normalized','outerposition',[0.05 0.1 0.9 0.85]);

set(gcf,'defaultAxesColorOrder',rbMap(NBK/2));
subplot(2,3,1); hold on;
plot(1e3*trueTime, 1e6*trueSol, 'linewidth',3,'color','k');
plot(1e3*trueTime, 1e6*squeeze(trajectories(:,2:2:end,9)),'linewidth',0.75);
xlim([1 10]);
xlabel('Time (ms)');
ylabel('[Ca^{2+}] (µM)');
title('Different BK @best Kappa');
cb = colorbar;
ylabel(cb,'BK');
colormap(gca,rbMap(NBK/2));
caxis([bk(1) bk(end)]);
set(gca,'fontsize',fontSize);

set(gcf,'defaultAxesColorOrder',rbMap(NBK/2));
subplot(2,3,4); hold on;
plot(1e3*trueTime, 1e6*(squeeze(trajectories(:,2:2:end,9))-trueSol),'linewidth',0.75);
xlim([1 10]);
xlabel('Time (ms)');
ylabel('Error [Ca^{2+}] (µM)');
title('Error Different BK @best Kappa');
cb = colorbar;
ylabel(cb,'BK');
colormap(gca,rbMap(NBK/2));
caxis([bk(1) bk(end)]);
set(gca,'fontsize',fontSize);



% Cost Function
subplot(2,3,2);
imagesc(ks,bk,costs);
xlabel('KS');
ylabel('bk');
title('Cost as function of {bk,ks}');
cb = colorbar;
ylabel(cb,'Error');
colormap(gca,flipud(hot()));
caxis([min(costs(:)) max(costs(:))]);
set(gca,'colorscale','log');
set(gca,'fontsize',fontSize);

subplot(2,3,5);
plot(bk,costs(:,9),'linewidth',1.5,'color','k');
ylim([5e-18 1e-16]);
set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('BK');
ylabel('Cost');
title('Cost for all BK at best Kappa');
set(gca,'fontsize',fontSize);

% Trajectories - as function of ks
set(gcf,'defaultAxesColorOrder',rbMap(NKS));
subplot(2,3,3); hold on;
plot(1e3*trueTime, 1e6*squeeze(trajectories(:,end,:)),'linewidth',0.75);
plot(1e3*trueTime, 1e6*trueSol, 'linewidth',3,'color','k');
plot(1e3*prm3Time, 1e6*prm3Sol(:,1),'linewidth',3,'color',[0 0.7 0]);
xlim([1 10]);
xlabel('Time (ms)');
ylabel('[Ca^{2+}] (µM)');
title('Different Kappa @best BK');
cb = colorbar;
ylabel(cb,'Kappa');
colormap(gca,rbMap(NKS));
caxis([ks(1) ks(end)]);
set(gca,'fontsize',fontSize);

% Peak Values
pks = squeeze(max(trajectories,[],1));
truePk = max(trueSol);
subplot(2,3,6);
imagesc(ks,bk,pks-truePk);
colorbar;
colormap(gca,redblue);
caxis([-1 1]*max(abs(pks(:)-truePk(:))))
xlabel('KS');
ylabel('bk');
title('Error in Peak Value');
set(gca,'fontsize',fontSize);

% print(gcf,'-painters',fullfile(fpath,'reportBKvsKS_simplification'),'-dpng');


%% Get Best Values and Compare with Real Kon/Koff/Bt Model

% --- include 3param model from reparameterized data ---
[~,bestIdx] = min(costs(:));
[bkIdx,ksIdx] = ind2sub(size(costs),bestIdx);
koff = bk(bkIdx)/ks(ksIdx);

NChecks = 100;
kd = logspace(-6,-2,NChecks);
kon = koff./kd;
bt = bk(bkIdx)./kon;

tVec = 0:0.0001:0.01;
prm3Sols = zeros(length(tVec),NChecks);
prm3NoSat = zeros(length(tVec),NChecks);

tspan = [0, 0.01];
for nc = 1:NChecks
    kappa = bt(nc) / (systemPrms.rest + koff/kon(nc));
    u0 = [systemPrms.rest; systemPrms.rest*kappa];
    prms = [kon(nc), koff, bt(nc)];
    cSolution = ode23s(@(t,y) massActionDynamics(t,y,systemPrms, prms), tspan, u0, odeOptions);
    cEvalSol = deval(cSolution,tVec);
    prm3Sols(:,nc) = cEvalSol(1,:);
    
    cSolution = ode23s(@(t,y) massActionDynamics_noSat(t,y,systemPrms, prms), tspan, u0, odeOptions);
    cEvalSol = deval(cSolution,tVec);
    prm3NoSat(:,nc) = cEvalSol(1,:);
end

f = figure(10); 
rMap = cat(2, linspace(0.2,1,NChecks)', zeros(NChecks,2));
set(gcf,'defaultAxesColorOrder',rMap);
set(gcf,'units','normalized','outerposition',[4.609375e-02 4.812500e-01 4.507813e-01 3.527778e-01]);

clf;
subplot(2,3,1);
plot(tVec,prm3Sols)

subplot(2,3,2);
plot(tVec,prm3NoSat);

subplot(2,3,3);
plot(tVec,prm3Sols ./ repmat(bt,length(tVec),1));


systemPrms.amp = 30e-12;
for nc = 1:NChecks
    kappa = bt(nc) / (systemPrms.rest + koff/kon(nc));
    u0 = [systemPrms.rest; systemPrms.rest*kappa];
    prms = [kon(nc), koff, bt(nc)];
    cSolution = ode23s(@(t,y) massActionDynamics(t,y,systemPrms, prms), tspan, u0, odeOptions);
    cEvalSol = deval(cSolution,tVec);
    prm3Sols(:,nc) = cEvalSol(1,:);
    
    cSolution = ode23s(@(t,y) massActionDynamics_noSat(t,y,systemPrms, prms), tspan, u0, odeOptions);
    cEvalSol = deval(cSolution,tVec);
    prm3NoSat(:,nc) = cEvalSol(1,:);
end
subplot(2,3,4);
plot(tVec,prm3Sols)

subplot(2,3,5);
plot(tVec,prm3NoSat);

subplot(2,3,6);
plot(tVec,prm3Sols ./ repmat(bt,length(tVec),1));







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
    curr = restCurrent + ica(t,prm.amp)/(96485*prm.v*prm.ks);
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
    curr = restCurrent + ica(t,systemPrms.amp)/(96485*systemPrms.v);
    extrusion = systemPrms.beta * y(1);
    
    % Association/Dissociation with quasi-instant buffer
    association = y(1)*(prms(3)-y(2))*prms(1);
    dissociation = y(2)*prms(2);
    bufferExchange = association - dissociation;
    
    % Derivative
    dy(1,1) = curr - extrusion - bufferExchange;
    dy(2,1) = bufferExchange;
end

% Gives solution with mass-action buffer - reparameterized to bk
function dy = massActionDynamicsRPRM(t,y,systemPrms,prms)
    % Parameters: [v,ks,beta,rest,amp]
    % bk,kappa is optimization parameter
    
    % y(1) = cafree
    % y(2) = boundBuffer
    
    % Current/Extrusion Term
    restCurrent = systemPrms.beta * systemPrms.rest;
    curr = restCurrent + ica(t,systemPrms.amp)/(96485*systemPrms.v);
    extrusion = systemPrms.beta * y(1);
    
    % Association/Dissociation with quasi-instant buffer
    association = y(1)*prms(1);
    dissociation = y(2)*prms(1)/prms(2);
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
    curr = restCurrent + ica(t,systemPrms.amp)/(96485*systemPrms.v);
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








