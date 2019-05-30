

% dydt = bufferDynamics_011(t,y,ica,p);

%% CRAZY FUCKING GRID

dt = 5e-8; % seconds
ds = 1000; % downsample factor
T = 0.3;
tvec = 0:dt*ds:T;
NT = length(tvec);

% Grid parameters
kon = 5*10^8; % Kon just use 1
kd = [400e-9 800e-9 5e-6 40e-6 100e-6]; % Kds
beta = 20*1000./[5:5:20 50]; % (Assuming an endogenous binding ratio of 20...)
rest = [50e-9 100e-9 500e-9]; % Molar Calcium resting
kdString = {'400nM' '800nM' '5µM' '40µM' '100µM'};

NON = length(kon); 
NOF = length(kd); 
NB = length(beta); 
NR = length(rest); 

% system parameters
p.v = (4/3) * pi * 0.4^3 * 1e-15; % volume expressed as L 
p.kon = kon(1); % On-Rate 1/M-s
p.koff = koff(1);
p.bconc = 300e-6; 
p.beta = beta(1);
p.rest = rest(1);

kappa = @(p,ca) p.bconc / (ca + p.koff/p.kon); % Inline for equilibrium binding ratio

% 10pA current for 2ms -- assume instantaneous 20 Kb endogenous buffer!!!
ica = @(t) (1/20)*2e-12 * exp(-(t-2e-3).^2./(0.55e-3).^2) .* (t >= 1e-3) .* (t <= 3e-3); 

makeGrid = true;
if makeGrid
    y = pacell([NON NOF NB NR],[NT 2],'nan');
    kappas = zeros(NON,NOF,NB,NR);
    dynamicKappa = pacell([NON NOF NB NR],[NT 1],'nan');
    for non = 1:NON
        p.kon = kon(non);
        for nof = 1:NOF
            p.koff = kd(nof) * kon(non);
            for nb =1:NB
                p.beta = beta(nb);
                for nr = 1:NR
                    fprintf(...
                        'K_{on}: %d/%d | K_{off}: %d/%d | Beta: %d/%d | Rest: %d/%d\n',non,NON,nof,NOF,nb,NB,nr,NR);

                    p.rest = rest(nr);
                    kappas(non,nof,nb,nr) = kappa(p,p.rest);
                    iState = [p.rest p.rest*kappa(p,p.rest)];
                    [~,y{non,nof,nb,nr}] = eulerapp(@(t,y) bufferDynamics_011(t,y,ica,p),[0 T],iState,dt,ds);
                    dynamicKappa{non,nof,nb,nr} = y{non,nof,nb,nr}(:,2) ./ y{non,nof,nb,nr}(:,1); % dynamic kappa
                end
            end
        end
    end
end

save('dynamicKappa_y','y','dynamicKappa','kappas','kon','kd','beta','rest','tvec')


  


%% -- some plotting --
xx = cell2mat(squeeze(dynamicKappa(1,:,2,2)));
yy = cell2mat(squeeze(cellfun(@(c) c(:,1), y(1,:,2,2),'uni',0)));

f = callFigs(1);
set(gcf,'outerposition',[0 0 0.6781 0.9000]);
subplot(2,1,1); hold on;
plotvc(1000*tvec,yy,f);
xlim([0 10]);
ylabel('[Ca^{2+}]');
title('Calcium Evoked by 2ms 2pA I_{Ca}');
set(gca,'yticklabel',readableMolar(get(gca,'ytick')));
set(gca,'fontsize',16);


subplot(2,1,2); hold on;
plotvc(1000*tvec,xx./repmat(xx(1,:),size(xx,1),1),f)
xlim([0 10]);
ylabel('Relative \kappa');
title('Dynamic \kappa re: K_{off}');
legend(kdString{:},'location','northeast');
set(gca,'fontsize',16);



%%
g = callFigs(2);
subplot(2,1,1); hold on;
bar(max(yy,[],1));


%%

dyKappa = cell2mat(squeeze(dynamicKappa(1,:,1,1)));
justCalcium = cell2mat(squeeze(cellfun(@(c) c(:,1), y(1,:,2,1),'uni',0)));
justBuffer = cell2mat(squeeze(cellfun(@(c) c(:,2), y(1,:,2,1),'uni',0)));

eqRatio = zeros(size(dyKappa));
for ikd = 1:NOF
    eqRatio(:,ikd) = (p.bconc ./ (kd(ikd) + justCalcium(:,ikd)));
end

f = callFigs(1);
subplot(2,1,1); hold on;
plotvc(1000*tvec,dyKappa./repmat(dyKappa(1,:),size(dyKappa,1),1),f)
xlim([0 10]);
ylabel('Relative \kappa');
title('Dynamic \kappa re: K_{off}');
legend(kdString{:},'location','northeast');
set(gca,'fontsize',16);

subplot(2,1,2); hold on;
plotvc(1000*tvec,dyKappa./eqRatio,f);
xlim([0 10]);
ylabel('Ratio dynamic:equilibrium \kappa');
title('Dynamic \kappa Ratio');
legend(kdString{:},'location','northeast');
set(gca,'fontsize',16);


%% Derivative Binding Ratio

h = 10;
fpt = tvec(1+2*h:end-2*h);
dyKappa = cell2mat(squeeze(dynamicKappa(1,:,1,1)));
derKappa = cell2mat(squeeze(cellfun(@(c) fivePointDer(c(:,2),h)./fivePointDer(c(:,1),h),y(1,:,2,1),'uni',0)));
eqRatio = zeros(size(dyKappa));
for ikd = 1:NOF
    eqRatio(:,ikd) = (p.bconc ./ (kd(ikd) + justCalcium(:,ikd)));
end

f = callFigs(1);
subplot(2,1,1); hold on;
plotvc(1000*tvec,dyKappa./repmat(dyKappa(1,:),size(dyKappa,1),1),f);
xlim([0 10]);
ylabel('d[BCa]/d[Ca]');
title('Derivative \kappa');
legend(kdString{:},'location','northeast');
set(gca,'fontsize',16);

subplot(2,1,2); hold on;
plotvc(1000*fpt,derKappa,f);
xlim([0 10]);
ylabel('Ratio dynamic:equilibrium \kappa');
title('Dynamic \kappa Ratio');
legend(kdString{:},'location','northeast');
set(gca,'fontsize',16);




%% -- plot decays so you can see the time constant --

xx = cell2mat(squeeze(dynamicKappa(1,:,2,2)));
yy = cell2mat(squeeze(cellfun(@(c) c(:,1), y(1,:,2,2),'uni',0)));
zz = cell2mat(squeeze(cellfun(@(c) c(:,2), y(1,:,2,2),'uni',0)));


yyBase = yy - repmat(yy(1,:),size(yy,1),1); 
yyNorm = yyBase./repmat(max(yyBase,[],1),size(yy,1),1);

zzBase = zz - repmat(zz(1,:),size(zz,1),1);
zzNorm = zzBase./repmat(max(zzBase,[],1),size(zz,1),1);

f = callFigs(1);
set(gcf,'outerposition',[0 0 0.6781 0.9000]);
subplot(2,1,1); hold on;
plotvc(1000*tvec,yyNorm,f);
xlim([0 200]);
ylabel('normalized [Ca^{2+}]');
title('Calcium Evoked by 2ms 2pA I_{Ca}');
legend(kdString,'location','northeast');
set(gca,'fontsize',16);

subplot(2,1,2); hold on;
plotvc(1000*tvec,zzNorm,f)
xlim([0 200]);
ylabel('Relative \kappa');
title('Dynamic \kappa re: K_{off}');
legend(kdString{:},'location','northeast');
set(gca,'fontsize',16);






