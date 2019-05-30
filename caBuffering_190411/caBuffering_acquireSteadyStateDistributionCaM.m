

uiwait();

p.konTC = 8.4 * 10^7;
p.konRC = 2.5 * 10^7;
p.koffTC = 2.6 * 10^3;
p.koffRC = 6.5;
p.konTN = 7.7 * 10^8;
p.konRN = 3.2 * 10^10;
p.koffTN = 1.6 * 10^5;
p.koffRN = 2.2 * 10^4;
p.legend = {'N_0C_0','N_1C_0','N_2C_0','N_0C_1','N_1C_1','N_2C_1','N_0C_2','N_1C_2','N_2C_2'};


initState = [50e-9; 150e-6; zeros(8,1)];
dt = [1e-6*ones(1,9) 1e-8*ones(1,6)];
ds = [1e4*ones(1,9) 1e6*ones(1,6)];
T = [5*ones(1,7), 0.5*ones(1,4), 0.05*ones(1,4)];

iCalcium = logspace(log10(10e-9),log10(0.001),15);

NIC = length(iCalcium);
state = zeros(NIC,9);
for ic = 1:NIC
    fprintf(1,'Calcium concentration %d/%d\n',ic,NIC);
    initState(1) = iCalcium(ic);
    [t,y,d] = eulerapp(@(t,y) calmodulinModel_HoldCalcium(t,y,p),[0 T(ic)],initState,dt(ic),ds(ic));
    state(ic,:) = y(end,2:end);
end


ss = struct();
ss.iCalcium = iCalcium;
ss.state = state;
ss.dt = dt;
ss.ds = ds;
ss.T = T;
ss.initStateCaM = initState(2:end);

% save('steadyStateCaM','ss');

%% -
recoverVars = false;
if recoverVars
   unpackStruct('steadyStateCaM');
end
    
%% Simple Plot 

statePlot = state./repmat(max(state,[],1),size(state,1),1);
axisPlot = max(state,[],1)./max(state(:));

subplot(5,1,[1 4]);
imagesc(state);
colormap('hot');
set(gca,'xticklabel',p.legend);

subplot(5,1,5);
imagesc(axisPlot);
colormap('hot');




%% Determine change in fluorescence from different starting points

hpath = '/Users/LandauLand/Documents/Research/SabatiniLab/modeling/caBuffering_190411';
load(fullfile(hpath,'steadyStateCaM'));
NIC = length(ss.iCalcium);

p.konTC = 8.4 * 10^7;
p.konRC = 2.5 * 10^7;
p.koffTC = 2.6 * 10^3;
p.koffRC = 6.5;
p.konTN = 7.7 * 10^8;
p.konRN = 3.2 * 10^10;
p.koffTN = 1.6 * 10^5;
p.koffRN = 2.2 * 10^4;
p.legend = {'N_0C_0','N_1C_0','N_2C_0','N_0C_1','N_1C_1','N_2C_1','N_0C_2','N_1C_2','N_2C_2'};
p.bound = [0 1 2 1 2 3 2 3 4];

caStep = 1e-6*[0.1 0.25 0.5 1 2 5 10];
NCS = length(caStep);

dt = 1e-6;
ds = 100;
T = 2;

response = cell(NIC,NCS);
for ic = 1:NIC
    for cs = 1:NCS
        fprintf('| IC: %d/%d | CS: %d/%d...\n',ic,NIC,cs,NCS);
        initState = [caStep(cs); ss.state(ic,:)'];
        [t,y,d] = eulerapp(@(t,y) calmodulinModel_HoldCalcium(t,y,p),[0 T],initState,dt,ds);
        response{ic,cs} = y;
    end
end


caStepResponse = struct();
caStepResponse.dt = dt;
caStepResponse.ds = ds;
caStepResponse.T = T;
caStepResponse.p = p;
caStepResponse.caStep = caStep;
caStepResponse.iCalcium = iCalcium;
caStepResponse.iState = ss.state;
caStepResponse.response = response;
% save('stepCalciumCaM','caStepResponse');


%% --- plots ---

ic = 6;
sc = 5;

f = figure(1);
subplot(1,2,1);
plotvc(t,cs.response{ic,sc}(:,3:end),f);
xlim([0 T]);
legend(p.legend(2:end),'location','southeast')


subplot(1,2,2);
bar(mean(cs.response{ic,sc}(round(end/2):end,3:end),1));
set(gca,'xticklabel',p.legend(2:end))


%% -- plot Binding --

totalBound = size(iCalcium);
for i = 1:length(iCalcium)
    totalBound(i) = sum(state(i,:).*p.bound) / (150e-6 * 4);
end
plot(iCalcium, totalBound);
line(xlim,[0.5 0.5]);
set(gca,'xticklabel',readableMolar(get(gca,'xtick')))



