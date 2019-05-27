

hpath = '/Users/LandauLand/Documents/Research/SabatiniLab/modeling/caBuffering_190411';

% General note on units:
% using standard units


%% Basic Model

dt = 1e-6; % seconds

% system parameters
p.v = (4/3) * pi * 0.4^3 * 1e-15; % volume expressed as L 
p.kon = 5*10^8; % On-Rate 1/M-s
p.koff = 800e-9 * p.kon; % Off-Rate 1/s
p.bconc = 0.0003; 
p.beta = 1/0.012; % 12ms time constant
p.rest = 5e-8; % Molar Calcium resting

kappa = p.bconc / (p.rest + p.koff/p.kon);

T = 1;
NT = T/dt;
tvec = dt:dt:T;

y = zeros(NT,2);
y(1,1) = p.rest;
y(1,2) = y(1,1) * kappa;

ica = zeros(NT,1);
iCenter = 0.102/dt;
iWidth = 0.001/dt;
iAmplitude = 10e-15;
ica(iCenter-iWidth*5:iCenter+iWidth*5) = iAmplitude * exp(-(-iWidth*5:iWidth*5).^2/(2*iWidth^2));



for t = 2:NT
    dydt = bufferDynamics_011(t,y(t-1,:),ica(t),p);
    y(t,:) = y(t-1,:) + dydt * dt;
end

dca = (y(:,1) - p.rest)/p.rest;
dff = (y(:,2) - p.rest*kappa)/(p.rest*kappa);

figure(1); hold on;
plot(tvec,dca,'k');
plot(tvec,dff,'r');
% xlim([tvec(NT/100) tvec(end)]);
legend('Calcium','Buffer','location','northeast');





%% Check compare kon and kd

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

konArray = 1*10.^(3:9);
kdArray = 1*10.^(-9:-4);

nKN = length(konArray);
nKD = length(kdArray);

msg = '';
fprintf(1,'\n\n\n| Running Model Comparing Kon and Kd...\n');

y = cell(nKN,nKD);
for kn = 1:nKN
    for kd = 1:nKD
        msg = cycleMessage(msg,sprintf('| Kon: %d/%d | Kd: %d/%d | Time Step: 0/100...',kn,nKN,kd,nKD));
        
        p.kon = konArray(kn);
        p.koff = kdArray(kd) * p.kon;
        kappa = p.bconc / (p.rest + p.koff/p.kon);
        
        y{kn,kd}(1,1) = p.rest*2;
        y{kn,kd}(1,2) = y{kn,kd}(1,1) * kappa;
        
        tCounter = true(1,10);
        for t = 2:NT
            if 10*t/NT >= find(tCounter,1,'first')
                percent = find(tCounter,1,'first')*10;
                msg = cycleMessage(msg,sprintf('| Kon: %d/%d | Kd: %d/%d | Time Step: %d/100...',kn,nKN,kd,nKD,percent));
                tCounter(find(tCounter,1,'first')) = false;
            end
            
            dydt = bufferDynamics_011(t,y{kn,kd}(t-1,:),0,p);
            y{kn,kd}(t,:) = y{kn,kd}(t-1,:) + dydt * dt;
        end
    end
end
fprintf(1,'| finished... \n');



%% Check compare kon and kd with ode45

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

konArray = 1*10.^(3:9);
kdArray = 1*10.^(-7:-4);

% konArray = 5*10^8;
% kdArray = 1*10^-5;

nKN = length(konArray);
nKD = length(kdArray);


options = odeset('NonNegative',[1 1],'MaxStep',dt*10);
fprintf(1,'\n\n| Running Model Comparing Kon and Kd...\n');
t = cell(nKN,nKD);
y = cell(nKN,nKD);
for kn = 1:nKN
    for kd = 1:nKD
        fprintf(1,'| Kon: %d/%d | Kd: %d/%d | working...\n',kn,nKN,kd,nKD);
        
        p.kon = konArray(kn);
        p.koff = kdArray(kd) * p.kon;
        kappa = p.bconc / (p.rest + p.koff/p.kon);
        
        initState = [p.rest p.rest*kappa(p.rest,p)]; 
        
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








































