

%% ---------------- now with same rate but different kernel --------------

T = 1e4;

trials = 100;

filter = logspace(log10(10),log10(1000),15);
rate = logspace(log10(1),log10(100),20);

NF = length(filter);
NR = length(rate);


sc = zeros(NR,NF);
msg = '';
for nf = 1:NF
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf('filter %d/%d...\n',nf,NF);
    fprintf(1,msg);
    fKernel = exp(-(0:T-1)/filter(nf));
    for nr = 1:NR
        currTrialCoeff = zeros(1,trials);
        for tr = 1:trials
            x1 = double(rand(T,1)<rate(nr)/1000);
            x2 = double(rand(T,1)<rate(nr)/1000);
            s1 = conv(x1,fKernel,'same');
            s2 = conv(x2,fKernel,'same');
            currTrialCoeff(tr) = corr(s1,s2);
        end
        sc(nr,nf) = mean(currTrialCoeff);
    end
end

% do special at 100ms filter
NT = 1000;
specialCoeff = zeros(NR,NT);
fKernel = exp(-(0:T-1)/100);
for nr = 1:NR
    for tr = 1:NT
        x1 = double(rand(T,1)<rate(nr)/1000);
        x2 = double(rand(T,1)<rate(nr)/1000);
        s1 = conv(x1,fKernel,'same');
        s2 = conv(x2,fKernel,'same');
        specialCoeff(nr,tr) = corr(s1,s2);
    end
end


%%
figure(1); clf;
set(gcf,'units','normalized','outerposition',[0.05 0.51 0.75 0.47]);

subplot(1,3,1);
imagesc(1:NF,1:NR,sc);
set(gca,'xtick',1:NF);
set(gca,'ytick',1:NR);
set(gca,'xticklabel',round(filter));
set(gca,'yticklabel',round(rate));
xlabel('time constant (ms)');
ylabel('rate (1/s)');
colorbar
set(gca,'fontsize',16);

cmap = [linspace(0.1,1,NR)',zeros(NR,2)];
subplot(1,3,2); hold on;
for nr = 1:NR
    plot(filter,sc(nr,:),'color',cmap(nr,:),'linewidth',1.5);
end
xlabel('time constant (ms)');
ylabel('c.coeff');
legend(cellfun(@(c) sprintf('%3d/s',c), num2cell(round(rate)),'uni',0),'location','southeast');
set(gca,'fontsize',16);

subplot(1,3,3); 
shadedErrorBar(rate,mean(specialCoeff,2),std(specialCoeff,[],2),{'color','k','linewidth',1.5});
xlim([1 100]);
ylim([0 1]);
xlabel('rate (1/s)');
ylabel('c.coeff');
title('w/ 100ms filter (+/- std)');
set(gca,'fontsize',16);



%% can be different rates but always 100ms

% do special at 100ms filter
NT = 100;
NR = 20;
rate = logspace(log10(1),log10(100),NR);
fKernel = exp(-(0:T-1)'/100);

rateVary100 = zeros(NR,NR,NT);
for nr1 = 1:NR
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf('rate %d/%d...\n',nr1,NR);
    fprintf(1,msg);
    
    for nr2 = 1:nr1
        for tr = 1:NT
            x1 = double(rand(T,1)<rate(nr1)/1000);
            x2 = double(rand(T,1)<rate(nr2)/1000);
            s1 = conv(x1,fKernel);
            s2 = conv(x2,fKernel);
            s1 = s1(1:T);
            s2 = s2(1:T);
            rateVary100(nr1,nr2,tr) = corr(s1,s2);
        end
    end
end

%%
figure(3); clf;
imagesc(1:NR,1:NR,mean(rateVary100,3));
set(gca,'xtick',1:NR);
set(gca,'ytick',1:NR);
set(gca,'xticklabel',round(10*rate)/10);
set(gca,'yticklabel',round(10*rate)/10);
colorbar;
caxis([0 1]);



%%
% T = 1e4;
% 
% 
% filter = logspace(log10(10),log10(1000),25);
% rate = logspace(log10(1),log10(100),10);
% 
% NF = length(filter);
% NR = length(rate);
% 
% 
% sc = zeros(NR,NR,NF);
% msg = '';
% for nr1 = 1:NR
%     fprintf(1,repmat('\b',1,length(msg)));
%     msg = sprintf('rate1 %d/%d...\n',nr1,NR);
%     fprintf(1,msg);
%     for nr2 = 1:nr1
%         x1 = double(rand(T,1)<rate(nr1)/1000);
%         x2 = double(rand(T,1)<rate(nr2)/1000);
%         for nf1 = 1:NF
%             fKernel = exp(-(0:T-1)/filter(nf1));
%             s1 = conv(x1,fKernel,'same');
%             s2 = conv(x2,fKernel,'same');
%             sc(nr1,nr2,nf1) = corr(s1,s2);
%         end
%     end
% end
% 
% 
% %%
% figure(1); clf;
% set(gcf,'units','normalized','outerposition',[0.05 0.55 0.64 0.37]);
% 
% subplot(1,3,1);
% imagesc(rate,rate,sc(:,:,1))
% title('Correlation w/ Filter 10ms');
% set(gca,'fontsize',16);
% 
% subplot(1,3,2);
% imagesc(rate,rate,sc(:,:,13));
% title('Correlation w/ Filter 100ms');
% set(gca,'fontsize',16);
% 
% subplot(1,3,3);
% imagesc(rate,rate,sc(:,:,25));
% title('Correlation w/ Filter 1000ms');
% set(gca,'fontsize',16);
% 







