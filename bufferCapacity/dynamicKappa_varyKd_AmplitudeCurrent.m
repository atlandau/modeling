

% dydt = bufferDynamics_011(t,y,ica,p);

%% CRAZY FUCKING GRID

dt = 5e-8; % seconds
ds = 1000; % downsample factor
T = 0.3;
tvec = 0:dt*ds:T;
NT = length(tvec);

% Grid parameters
kon = 5*10^8; % Kon just use 1
kd = 5e-6; %[400e-9 800e-9 5e-6 40e-6 100e-6]; % Kds
beta = 20*1000/12; % (Assuming an endogenous binding ratio of 20...)
rest = 75e-9; % Molar Calcium resting
kdString = {'400nM' '800nM' '5然' '40然' '100然'};
amplitude = 5e-12; %[0.5e-12 2e-12 5e-12 20e-12 50e-12];

NK = length(kd); 
NA = length(amplitude);

% system parameters
p.v = (4/3) * pi * 0.4^3 * 1e-15; % volume expressed as L 
p.kon = kon; % On-Rate 1/M-s
p.koff = [];
p.bconc = 300e-6; 
p.beta = beta;
p.rest = rest;
p.amplitude = amplitude(1); 

kappa = @(p,ca) p.bconc / (ca + p.koff/p.kon); % Inline for equilibrium binding ratio

% 10pA current for 2ms -- assume instantaneous 20 Kb endogenous buffer!!!
ica = @(t,amp) (1/20) * amp * exp(-(t-2e-3).^2./(0.55e-3).^2) .* (t >= 1e-3) .* (t <= 3e-3); 

makeGrid = true;
if makeGrid
    y = pacell([NK NA],[NT 2],'nan');
    kappas = zeros([NK NA]);
    dynamicKappa = pacell([NK NA],[NT 1],'nan');
    for nk = 1:NK
        p.koff = p.kon * kd(nk);
        for na = 1:NA
            p.amplitude = amplitude(na);
                fprintf(...
                    'Kd: %d/%d | Amplitude: %d/%d\n',nk,NK,na,NA);

            kappas(nk,na) = kappa(p,p.rest);
            iState = [p.rest p.rest*kappa(p,p.rest)];
            [~,y{nk,na}] = euclidapp(@(t,y) bufferDynamics_011(t,y,ica,p),[0 T],iState,dt,ds);
            dynamicKappa{nk,na} = y{nk,na}(:,2) ./ y{nk,na}(:,1); % dynamic kappa
        end
    end
end

% save('dynamicKappa_Kd_Amplitude','y','dynamicKappa','kappas','kon','kd','beta','rest','tvec')



%% -- analysis --

minDynKappa = cellfun(@min, dynamicKappa, 'uni', 1);
bufferCapacity = @(bt,kd,ca0,caAP) kd*bt / (kd + ca0) ./ (kd + caAP);

ca0 = rest;
caAP = cellfun(@(c) max(c(:,1))-ca0,y,'uni',1);
estimateKappa = zeros(NK,NA,2);
for nk = 1:NK
    for na = 1:NA
        deltaCaAP = max(y{nk,na}(:,1)) - ca0;
        estimateKappa(nk,na,1) = kd(nk) * p.bconc / (kd(nk) + ca0) / (kd(nk) + deltaCaAP);
        estimateKappa(nk,na,2) = kd(nk) * p.bconc / (kd(nk) + ca0) / (kd(nk) + ca0);
    end
end

callFigs(11);
subplot(2,3,1); hold on;
imagesc(estimateKappa(:,:,1)./kappas);
set(gca,'ydir','normal');
colormap('hot');
colorbar;
xlim([0.5 NA+0.5]);
ylim([0.5 NK+0.5]);
set(gca,'xtick',1:NA); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NK);
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
xlim([0.5 NA+0.5]);
ylim([0.5 NK+0.5]);
set(gca,'xtick',1:NA); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NK);
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
xlim([0.5 NA+0.5]);
ylim([0.5 NK+0.5]);
set(gca,'xtick',1:NA); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NK);
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
xlim([0.5 NA+0.5]);
ylim([0.5 NK+0.5]);
set(gca,'xtick',1:NA); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NK);
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
xlim([0.5 NA+0.5]);
ylim([0.5 NK+0.5]);
set(gca,'xtick',1:NA); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NK);
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
xlim([0.5 NA+0.5]);
ylim([0.5 NK+0.5]);
set(gca,'xtick',1:NA); 
set(gca,'xticklabel',round(amplitude*1e12*10)/10);
set(gca,'ytick',1:NK);
set(gca,'yticklabel',readableMolar(kd));
xlabel('Amplitude (pA)');
ylabel('Kd');
ylabel(cb,'Resting \kappa');
title('Resting \kappa');
set(gca,'fontsize',16);

%% -- show calcium and buffer dynamics, raw and normalized, vary amplitude calcium influx --

xLimConcentrations = [0 30];
xLimKappas = [0 30];

kdIDX = kd==5e-6; % Grab all amplitudes from 5然 Kd
ca = cell2mat(cellfun(@(c) c(:,1), y(kdIDX,:), 'uni', 0));  % Calcium 
bf = cell2mat(cellfun(@(c) c(:,2), y(kdIDX,:), 'uni', 0));  % Buffer
dk = cell2mat(dynamicKappa(kdIDX,:)); % Dynamic kappa

kappaKD = @(bt,ca,kd) bt ./ (ca + kd);
eqDK = cellfun(@(c,kd) kappaKD(p.bconc, c(:,1), kd), y, num2cell(repmat(kd(:),1,NA)),'uni',0);
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
legend(cellfun(@(c) sprintf('%.1f pA',c*1e12),num2cell(amplitude),'uni',0),'location','northeast');
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
plotvc(1000*tvec, cell2mat(eqDK(kdIDX,:)), f, {'''linewidth''','0.5'});
xlim(xLimKappas);
ylabel('Buffer Ratio');
title('Dynamic vs. Equilibrium \kappa');
set(gca,'fontsize',16);

subplot(2,3,6); hold on;
plotvc(1000*tvec, cell2mat(dkRatio(kdIDX,:)), f);
xlim(xLimKappas);
ylabel('Ratio');
title('Dynamic/Equilibrium \kappa');
set(gca,'fontsize',16);































