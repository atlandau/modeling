

hpath = '/Users/LandauLand/Documents/Research/SabatiniLab/modeling/caBuffering_190411';


%% Compare Kon Kd / ode45 / stepResponse of calcium

dt = 1e-6; 

T = 1; 
NT = T/dt; 
tvec = dt:dt:T; 

% system parameters
p.v = (4/3) * pi * 0.4^3 * 1e-15; % volume expressed as L 
p.kon = 1*10^8; % On-Rate 1/M-s
p.koff = 0.00005 * p.kon; % Off-Rate 1/s
p.bconc = 0.0003; % 0 M
p.beta = 1/0.012; % 12ms time constant
p.rest = 5e-8; % Molar Calcium resting

kappa = @(ca,p) p.bconc / (ca + p.koff/p.kon); % inline function defining equilibrium binding ratio at a given [Ca]

konArray = 1*10.^(5:9);
kdArray = 1*10.^(-7:-4);

nKN = length(konArray);
nKD = length(kdArray);


options = odeset('NonNegative',[1 1],'MaxStep',dt*1);
fprintf(1,'\n\n| Running Model Comparing Kon and Kd...\n');
t = cell(nKN,nKD);
y = cell(nKN,nKD);
for kn = 1:nKN
    for kd = 1:nKD
        fprintf(1,'| Kon: %d/%d | Kd: %d/%d | working...\n',kn,nKN,kd,nKD);
        
        p.kon = konArray(kn);
        p.koff = kdArray(kd) * p.kon;
        
        initCa = p.rest*10;
        initState = [initCa initCa*kappa(initCa,p)]; 
        
        ica = zeros(NT,1);
        iCenter = 0.052/dt;
        iWidth = 0.001/dt;
        iAmplitude = 10e-15;
        ica(iCenter-iWidth*5:iCenter+iWidth*5) = iAmplitude * exp(-(-iWidth*5:iWidth*5).^2/(2*iWidth^2));
        
        currentVector = [tvec(:) ica(:)];
        [ct,cy] = ode45(@(t,y) bufferDynamics_011(t,y,currentVector,p),[0 T],initState,options);
        
        % Save into cell array
        t{kn,kd} = ct;
        y{kn,kd} = [(cy(:,1) - p.rest)/p.rest, (cy(:,2) - p.rest*kappa(p.rest,p))/(p.rest*kappa(p.rest,p))];
    end
end
fprintf(1,'| finished... \n\n');





%% --

% buffDynamics.p = p;
% buffDynamics.konArray = konArray;
% buffDynamics.kdArray = kdArray;
% buffDynamics.options = options;
% buffDynamics.kappa = kappa;
% buffDynamics.ica = currentVector;
% buffDynamics.t = t;
% buffDynamics.y = y;
% 
% save('bufferDynamicsResults_190412','buffDynamics');


ts = find(tvec>=0.1,1);







%%

dt = 1e-8; 

T = dt*3000;
NT = round(T/dt);
T = dt*NT;
tvec = dt:dt:T;

% system parameters
p.v = (4/3) * pi * 0.4^3 * 1e-15; % volume expressed as L 
p.kon = 5*10^8; % On-Rate 1/M-s
p.koff = 0.0000008 * p.kon; % Off-Rate 1/s
p.bconc = 0.0003; % 0 M
p.beta = 0;%1000/12; % 12ms time constant
p.rest = 5e-8; % Molar Calcium resting

% kappa = @(ca,p) p.kon*p.bconc / (p.kon*ca - p.koff);
kappa = @(ca,p) p.bconc / (ca + p.koff/p.kon); % inline function defining equilibrium binding ratio at a given [Ca]

options = odeset('NonNegative',[1 1],'MaxStep',dt);

initCa = p.rest*10;
initState = [initCa initCa*kappa(initCa,p)]; 
        
iCenter = 0;
iWidth = dt*2;
iAmplitude = 10e-15;

ica = @(t) iAmplitude * exp(-(t-iCenter).^2/(2*iWidth^2)); % inline to define calcium current
icaStep = @(t) initCa*2*p.v*96485 * (t>iCenter-dt/2 & t<iCenter+dt/2);

[ct,cy] = ode45(@(t,y) bufferDynamics_011(t,y,icaStep,p),[0 T],initState,options);
dy = [(cy(:,1) - initCa)/initCa, (cy(:,2) - initCa*kappa(initCa,p))/(initCa*kappa(initCa,p))];


figure(1); clf;
subplot(2,1,1);
plot(cy(:,1));
subplot(2,1,2);
plot(cy(:,2));



%%

[tt,yy,dd] = eulerapp(@(t,y) bufferDynamics_011(t,y,icaStep,p),[0 T],initState,dt/10);

toPlot = yy;

stIdx = 1;

figure(1); clf;
subplot(2,1,1);
plot(tt(stIdx:end)*1000000,yy(stIdx:end,1));  
subplot(2,1,2);
plot(tt(stIdx:end)*1000000,yy(stIdx:end,2));




%% -----
p.konTC = 7.9 * 10^7;
p.konRC = 7.4 * 10^7;
p.koffTC = 3.4 * 10^3;
p.koffRC = 1.2;
p.konTN = 8.9 * 10^8;
p.konRN = 10.5 * 10^10;
p.koffTN = 5.2 * 10^5;
p.koffRN = 4.3 * 10^4;
p.legend = {'N_0C_0','N_1C_0','N_2C_0','N_0C_1','N_1C_1','N_2C_1','N_0C_2','N_1C_2','N_2C_2'};

initState = [50e-9; 150e-6; zeros(8,1)];

dt = 1e-6;
T = 5;
[t,y,d] = eulerapp(@(t,y) calmodulinModel_HoldCalcium(t,y,p),[0 T],initState,dt);


f = figure;
subplot(1,2,1);
plotvc(t,y(:,3:end),f);
xlim([0 5]);
legend(p.legend(2:end),'location','southeast')
subplot(1,2,2);
bar(mean(y(end/2:end,3:end),1));
set(gca,'xticklabel',p.legend(2:end))




%     y(2) = N0_C0
%     y(3) = N1_C0
%     y(4) = N2_C0
%     y(5) = N0_C1
%     y(6) = N1_C1
%     y(7) = N2_C1
%     y(8) = N0_C2
%     y(9) = N1_C2
%     y(10) = N2_C2






%% -----
p.konTC = 7.9 * 10^7;
p.konRC = 7.4 * 10^7;
p.koffTC = 3.4 * 10^3;
p.koffRC = 1.2;
p.konTN = 8.9 * 10^8;
p.konRN = 10.5 * 10^10;
p.koffTN = 5.2 * 10^5;
p.koffRN = 4.3 * 10^4;
p.legend = {'N_0C_0','N_1C_0','N_2C_0','N_0C_1','N_1C_1','N_2C_1','N_0C_2','N_1C_2','N_2C_2'};

initState = [iCalcium(12); 150e-6; zeros(8,1)];

dt = 1e-8;
ds = 100000;
T = 0.01;
[t,y,d] = eulerapp(@(t,y) calmodulinModel_HoldCalcium(t,y,p),[0 T],initState,dt,ds);


f = figure;
subplot(1,2,1);
plotvc(t,y(:,3:end),f);
xlim([0 T]);
legend(p.legend(2:end),'location','southeast')
subplot(1,2,2);
bar(mean(y(round(end/2):end,3:end),1));
set(gca,'xticklabel',p.legend(2:end))







