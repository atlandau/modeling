function out = icaTrain(t,amp,onsetTimes)
    % function giving calcium current with fixed temporal properties and
    % user-defined amplitude
    
    nInputs = length(onsetTimes);
    out = zeros(size(t));
    for n = 1:nInputs
        relativeTime = t - onsetTimes(n);
        
        current = amp*exp(-(relativeTime-1e-3).^2 / (0.55e-3)^2);
        offset = amp*exp(-(1e-3).^2 / (0.55e-3)^2);
        adjustedCurrent = (current - offset) .* (relativeTime>=0).*(relativeTime<=2e-3);
        out = out + adjustedCurrent.*(adjustedCurrent>=0);
    end
end