

% Using the following model: dydt = bufferDynamics_011(t,y,ica,p);

% Time Parameters
dt = 1e-6; % seconds
ds = 100; % downsample factor
T = 1;
tvec = 0:dt*ds:T;
NT = length(tvec);


% Buffer Parameters
Ks = 20; %instantaneous endogenous kappa
gamma = 21/0.012; % rate of extrusion (1/second)
beta = gamma/Ks; % adjusted rate
rest = 75e-9; % resting calcium 
amplitude = 1e-12; % maximum pA calcium current

kon = 5*10^8; % on rate
kd = 2.3e-6; % kd (assumed)
koff = kon*kd;
bconc = 0;%40e-6; %40e-6:40e-6:300e-6; % concentration array

NC = length(bconc);

% system parameters
p.v = (4/3) * pi * 0.4^3 * 1e-15; % volume expressed as L 
p.ks = Ks;
p.kon = kon; % On-Rate 1/M-s
p.koff = koff;
p.bconc = bconc;
p.gamma = gamma;
p.rest = rest;
p.amplitude = amplitude; 

kappa = @(p,ca) p.bconc / (ca + p.koff/p.kon); % Inline for equilibrium binding ratio

% To estimate kappa inaccurately
kappaEst0 = @(kd,bt,ca0) kd*bt / (kd+ca0) / (kd+ca0);
kappaEstAP = @(kd,bt,ca0,caAP) kd*bt / (kd+ca0) / (kd+caAP);

% 10pA current for 2ms -- assume instantaneous 20 Kb endogenous buffer!!!
ica = @(t,amp,ks) (1/ks) * amp * exp(-(t-2e-3).^2./(0.55e-3).^2) .* (t >= 1e-3) .* (t <= 3e-3); 

startFit = find(tvec>=0.004,1);
endFit = find(tvec>=0.200,1);
fitidx = startFit:endFit;
fittime = tvec(fitidx)';

makeGrid = true;
if makeGrid
    y = pacell([1 NC],[NT 2],'nan');
    kappas = zeros(1,NC);
    dynamicKappa = pacell([1 NC],[NT 1],'nan');
    fits = cell(1,NC);
    tc = zeros(1,NC);
    kappaEstimate = zeros(2,NC,1);
    for nc = 1:NC
        for type=1:1
            fprintf('Type: %d | Concentration: %d/%d\n',type,nc,NC);
            
            p.bconc = bconc(nc); % only keep buffer if we're in type 1
            
            kappas(1,nc) = kappa(p,p.rest);
            iState = [p.rest p.rest*kappa(p,p.rest)];
            [~,y{type,nc}] = eulerapp(@(t,y) bufferDynamics_011(t,y,ica,p),[0 T],iState,dt,ds);
            
            dynamicKappa{type,nc} = y{1,nc}(:,2) ./ y{1,nc}(:,1); % dynamic kappa

            % Analysis
            kappaEstimate(1,nc,type) = kappaEst0(kd,p.bconc,p.rest); % estimate kappa
            kappaEstimate(2,nc,type) = kappaEstAP(kd,p.bconc,p.rest,max(y{1,nc}(:,1))); % estimate kappa

            d2fit = y{type,nc}(fitidx,2) - y{type,nc}(end,2);
            fits{type,nc} = fit(fittime,d2fit,'exp1'); % fits
            tc(nc) = -1000*1/fits{type,nc}.b; % time constant in ms
        end
    end
end

% save('kdvsamplitude','y','dynamicKappa','kappas','kon','kd','beta','rest','tvec')





%% -- show calcium and buffer dynamics, raw and normalized, vary amplitude calcium influx --

xLimConcentrations = [0 200];
xLimKappas = [0 200];

ca = cell2mat(cellfun(@(c) c(:,1), y, 'uni', 0));  % Calcium 
bf = cell2mat(cellfun(@(c) c(:,2), y, 'uni', 0));  % Buffer
dk = cell2mat(dynamicKappa); % Dynamic kappa

kappaKD = @(bt,ca,kd) bt ./ (ca + kd);
eqDK = cellfun(@(c,bt) kappaKD(bt, c(:,1), kd), y, num2cell(bconc),'uni',0);
dkRatio = cellfun(@(dk,eq) dk./eq, dynamicKappa, eqDK, 'uni', 0);

normCA = norman(ca);
normBF = norman(bf);
normDK = norman(dk);

f = callFigs(1);
subplot(2,3,1);
plotvc(1000*tvec,ca,f)
xlim(xLimConcentrations);
ylabel('Calcium [M]');
title('Calcium');
legend(cellfun(@(c) sprintf('%.1f µM',c*1e6),num2cell(bconc),'uni',0),'location','northeast');
set(gca,'fontsize',16);
% set(gca,'yticklabel',readableMolar(get(gca,'ytick')));

subplot(2,3,2);
plotvc(1000*tvec,normCA,f)
xlim(xLimConcentrations);
ylabel('Normalized Calcium');
title('Calcium');
set(gca,'fontsize',16);

subplot(2,3,4);
plotvc(1000*tvec,bf,f)
xlim(xLimConcentrations);
ylabel('Buffer [M]');
title('Buffer');
set(gca,'fontsize',16);
% set(gca,'yticklabel',readableMolar(get(gca,'ytick')));

subplot(2,3,5);
plotvc(1000*tvec,normBF,f)
xlim(xLimConcentrations);
ylabel('Normalized Buffer');
title('Buffer');
set(gca,'fontsize',16);

subplot(2,3,3); hold on;
plotvc(1000*tvec, dk, f, {'''linewidth''','2'});
plotvc(1000*tvec, cell2mat(eqDK), f, {'''linewidth''','0.5'});
xlim(xLimKappas);
ylabel('Buffer Ratio');
title('Dynamic vs. Equilibrium \kappa');
set(gca,'fontsize',16);

subplot(2,3,6); hold on;
plotvc(1000*tvec, cell2mat(dkRatio), f);
xlim(xLimKappas);
ylabel('Ratio');
title('Dynamic/Equilibrium \kappa');
set(gca,'fontsize',16);




%% -- do kappa estimation --

figure(1); clf;
hold on;
plot(kappaEstimate(1,:),tc,'color','k','marker','o','markerfacecolor','k','linestyle','none');
plot(kappaEstimate(2,:),tc,'color','r','marker','o','markerfacecolor','r','linestyle','none');

ke1 = [kappaEstimate(1,2:end)', ones(NC-1,1)] \ tc(2:end)';
ke2 = [kappaEstimate(2,2:end)', ones(NC-1,1)] \ tc(2:end)';

line([-40 200],ke1(1)*[-40 200]+ke1(2),'color','k');
line([-40 200],ke2(1)*[-40 200]+ke2(2),'color','r');

% xlim([-30 150]);
% ylim([0 100]);





%% -- analysis --

minDynKappa = cellfun(@min, dynamicKappa, 'uni', 1);
bufferCapacity = @(bt,kd,ca0,caAP) kd*bt / (kd + ca0) ./ (kd + caAP);

ca0 = rest;
caAP = cellfun(@(c) max(c(:,1))-ca0,y,'uni',1);
estimateKappa = zeros(NKD,NC,2);
for nkd = 1:NKD
    for nc = 1:NC
        deltaCaAP = max(y{nkd,nc}(:,1)) - ca0;
        estimateKappa(nkd,nc,1) = kd(nkd) * p.bconc / (kd(nkd) + ca0) / (kd(nkd) + deltaCaAP);
        estimateKappa(nkd,nc,2) = kd(nkd) * p.bconc / (kd(nkd) + ca0) / (kd(nkd) + ca0);
    end
end

callFigs(11);
subplot(2,3,1); hold on;
imagesc(estimateKappa(:,:,1)./kappas);
set(gca,'ydir','normal');
colormap('hot');
colorbar;
xlim([0.5 NC+0.5]);
ylim([0.5 NKD+0.5]);
set(gca,'xtick',1:NC); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NKD);
set(gca,'yticklabel',readableMolar(kd));
xlabel('Amplitude (pA)');
ylabel('Kd');
title('\kappa(ca0,caAP) / BaseKappa');
set(gca,'fontsize',16);

subplot(2,3,2); hold on;
imagesc(estimateKappa(:,:,2)./kappas);
set(gca,'ydir','normal');
colormap('hot');
colorbar;
xlim([0.5 NC+0.5]);
ylim([0.5 NKD+0.5]);
set(gca,'xtick',1:NC); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NKD);
set(gca,'yticklabel',readableMolar(kd));
xlabel('Amplitude (pA)');
ylabel('Kd');
title('\kappa(ca0,ca0) / BaseKappa');
set(gca,'fontsize',16);


subplot(2,3,4); hold on;
imagesc(estimateKappa(:,:,1)./minDynKappa);
set(gca,'ydir','normal');
colormap('hot');
colorbar;
xlim([0.5 NC+0.5]);
ylim([0.5 NKD+0.5]);
set(gca,'xtick',1:NC); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NKD);
set(gca,'yticklabel',readableMolar(kd));
xlabel('Amplitude (pA)');
ylabel('Kd');
title('\kappa(ca0,caAP) / minDynKappa');
set(gca,'fontsize',16);

subplot(2,3,5); hold on;
imagesc(estimateKappa(:,:,2)./minDynKappa);
set(gca,'ydir','normal');
colormap('hot');
colorbar;
xlim([0.5 NC+0.5]);
ylim([0.5 NKD+0.5]);
set(gca,'xtick',1:NC); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NKD);
set(gca,'yticklabel',readableMolar(kd));
xlabel('Amplitude (pA)');
ylabel('Kd');
title('\kappa(ca0,ca0) / minDynKappa');
set(gca,'fontsize',16);

subplot(2,3,3);
imagesc(cellfun(@(c) max(c(:,1)), y, 'uni', 1));
set(gca,'ydir','normal');
colormap('hot');
cb = colorbar;
xlim([0.5 NC+0.5]);
ylim([0.5 NKD+0.5]);
set(gca,'xtick',1:NC); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NKD);
set(gca,'yticklabel',readableMolar(kd));
xlabel('Amplitude (pA)');
ylabel('Kd');
ylabel(cb,'Max [Ca^{2+}]');
cb.TickLabels = readableMolar(cb.Ticks);
title('Max Calcium Concentration');
set(gca,'fontsize',16);

subplot(2,3,6);
imagesc(kappas);
set(gca,'ydir','normal');
colormap('hot');
cb = colorbar;
xlim([0.5 NC+0.5]);
ylim([0.5 NKD+0.5]);
set(gca,'xtick',1:NC); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NKD);
set(gca,'yticklabel',readableMolar(kd));
xlabel('Amplitude (pA)');
ylabel('Kd');
ylabel(cb,'Resting \kappa');
title('Resting \kappa');
set(gca,'fontsize',16);









