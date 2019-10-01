
tauNMDAR = 16; %ms
tauCa = 16; % ms

dt = 0.01; % ms
T = 150; % ms
tv = 0:dt:T; 
NT = length(tv); % number samples


nmda = exp(-tv/tauNMDAR);
calcium = exp(-tv/tauNMDAR);

ny = fft(nmda);
cy = fft(calcium);

bothy = ny.*cy;
both = real(ifft(bothy));

hold on; cla;
plot(tv, nmda,'k');
plot(tv,calcium,'r');
plot(tv,both./max(both),'b');









