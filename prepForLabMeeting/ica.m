function out = ica(t,amp)
    % Starts at 0ms

    % function giving calcium current with fixed temporal properties and
    % user-defined amplitude
    current = amp*exp(-(t-1e-3).^2 / (0.55e-3)^2);
    offset = amp*exp(-(1e-3).^2 / (0.55e-3)^2);
    adjustedCurrent = (current - offset) .* (t>=0).*(t<=2e-3);
    out = adjustedCurrent .* (adjustedCurrent>=0);
end