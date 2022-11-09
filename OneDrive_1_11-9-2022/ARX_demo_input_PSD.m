% ARX demo
%
% System identification using ARX model - input power spectrum

% input sequence
T = 16*4096;
time = 0:T-1;
p_switch = 0.1; % show low/high band properties of input signal
UU = ones(1,T);
for t=2:T,
    if rand(1)<p_switch, UU(t)=-UU(t-1); else UU(t)=UU(t-1); end
end
figure(1)
plot(0:500,UU(1:501))
figure(2)
[Ruu, lags]= xcov(UU,50);
stem(lags,Ruu)
figure(3)
[Puu, W]=pwelch(UU);
semilogx(W,log10(Puu))
grid on
hold on



