% low-frequency
% white noise
% high-frequency
function [UU] = createSignal(p_switch, n_samples)
    % input sequence - white noise
    UU = ones(n_samples,1);
    for t=2:(n_samples-1)
        if rand(1)<p_switch, UU(t)=-UU(t-1); else UU(t)=UU(t-1); end
    end
end