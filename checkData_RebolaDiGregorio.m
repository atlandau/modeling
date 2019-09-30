
% Rate Terms
diffusionLimited = 5e8;
egta_on = 1e7;
egta_off = 0.7;
ogb5n_on = diffusionLimited;
ogb5n_off = diffusionLimited * 20e-6; % Kon * Kd
fluo5f_on = diffusionLimited;
fluo5f_off = diffusionLimited * 2.3e-6; % Kon * Kd

% EGTA
egta.on = egta_on;
egta.off = egta_off;
egta.concentration = 0.1e-3;
% OGB5N
ogb.on = ogb5n_on;
ogb.off = ogb5n_off;
ogb.concentration = 100e-6;
% Fluo5f
fluo5f.on = fluo5f_on;
fluo5f.off = fluo5f_off;
fluo5f.concentration = 100e-6;

% Calcium
calciumConcentration = 20e-6; 

% Construct and Run Problem
tspan = [0,0.1];
tolerance = 1e-10;
odeOptions = odeset('RelTol',tolerance,'AbsTol',tolerance);
dt = 0.0001;
tvec = tspan(1):dt:tspan(2);

iState = [calciumConcentration 0 0]';
odeProblem = @(t,y) twoBufferEquilibrium(t,y,egta,ogb);
[t,y] = ode23s(odeProblem,tspan,iState,odeOptions);

cols = 'kbr';
figure(1); clf;
hold on;
for i = 1:3
    plot(t,y(:,i),'color',cols(i),'linewidth',1.5);
end

logf = floor(log10(y(end,1)));
fprintf(1,'Resting Calcium: %.2f*10^%d\n',y(end,1)/(10^logf),logf)


%% --- functions used for equilibrium conditions --
function dy = twoBufferEquilibrium(~,y,buffer1,buffer2)
    % d[Ca]/dt = 
    % -kon1*[ca]*[buffer1] 
    % -kon2*[ca]*[buffer2]
    % +koff1*[Ca-Buffer1] 
    % +koff2*[Ca-Buffer2];

    freeBuffer1 = buffer1.concentration - y(2);
    freeBuffer2 = buffer2.concentration - y(3);
    b1Exchange = buffer1.on*y(1)*freeBuffer1 - buffer1.off*y(2);
    b2Exchange = buffer2.on*y(1)*freeBuffer2 - buffer2.off*y(3);

    dy(1,1) = -(b1Exchange + b2Exchange);
    dy(2,1) = b1Exchange;
    dy(3,1) = b2Exchange;
end








