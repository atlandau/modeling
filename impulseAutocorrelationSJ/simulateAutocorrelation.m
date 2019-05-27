% https://stats.stackexchange.com/questions/114347/how-is-the-impulse-response-function-of-a-given-system-related-to-the-autocorrel

T = 100001;
x = double(rand(T,1)<0.01);
k = exp(-(1:T)'/100);
y = conv(x,k);
y = y(1:T);
y = y + 0*randn(T,1);

% Autocorrelations
Rx = xcorr(x,x); % autocorrelation of the input
Rh = xcorr(k,k); % autocorrelation of the impulse response
Ry = xcorr(y,y); % autocorrelation of the output

ffEstimateRy = fftshift(ifft(fft(Rh).*fft(Ry))); % Rx * Rh (convolution of Rx and Rh)
normFFRy = ffEstimateRy .* mean(Ry./ffEstimateRy); % Normalized(Rx * Ry) because they are scaled differently

figure;(1);
subplot(2,1,1);
hold on;
plot(x,'k'); % input signal
plot(y,'r'); % output signal

subplot(2,1,2);
hold on;
plot(normFFRy,'k'); % Estimate of Autocorrelation - Normalized (Rx * Ry)
plot(Ry,'r'); % Autocorrelation of Output
