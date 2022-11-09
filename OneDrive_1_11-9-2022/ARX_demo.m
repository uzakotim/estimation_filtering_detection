% ARX demo
%
% System identification using ARX model

% transfer function with unit gain
a = conv([1 -0.7],[1 -0.9]); % a = [1 -1.6 0.63]
b = 0.03*[3 2 1]/6;          % unit gain
%
% input sequence
T = 1001;
time = 0:T-1;
p_switch = 0.1;      % pseudorandom binary signal
UU = 1*ones(1,T);    % set amplitude - S/N ratio
EE = randn(1,T);
for t=2:T,
    if rand(1)<p_switch, UU(t)=-UU(t-1); else UU(t)=UU(t-1); end
end
YY =filter(b,a,UU);
YYM =YY + filter(0.3,a,EE) ; % filter gain 10
%
% plot data
figure(1)
subplot(3,1,1)
plot(time,UU)
ylabel('u(t)')
title('ARX model data')
subplot(3,1,2)
plot(time,EE)
ylabel('e(t)')
subplot(3,1,3)
plot(time,YY,time,YYM)
ylabel('y(t)')
xlabel('time')
%
% initialize ARX model statistics for LDFIL
theta = zeros(5,1);
d_theta = 2*ones(5,1);
L_theta = eye(5,5);
nu = 2;
nus2 = 1;
%
stheta = [theta, d_theta, L_theta'];
ssigma = [nu, nus2];
%
AA = zeros(2,T);
BB = zeros(3,T);
SS2 = zeros(1,T);
for t=3:T,
    % create regressor z'=[dy,d2y, u, du,d2u]
    data = [-YYM(t-1),-YYM(t-2),UU(t),UU(t-1),UU(t-2),YYM(t)];
    % run LDFIL step, forgetting phi=1
    [stheta, ssigma, k, eps, dy] = ldfil(stheta, ssigma, data, 1);
    % select estimated parameters
    theta = stheta(:,1);
    AA(:,t)=theta(1:2);
    BB(:,t)=theta(3:5);
    SS2(t) = ssigma(2)/ssigma(1);
end
figure(2)
subplot(3,1,1)
plot(time,AA)
hold on
plot([0,T-1],[a(2),a(2)],'b:')
plot([0,T-1],[a(3),a(3)],'r:')
hold off
ylabel('a(t)')
title('ARX model parameters')
subplot(3,1,2)
plot(time,BB)
hold on
plot([0,T-1],[b(1),b(1)],'b:')
plot([0,T-1],[b(2),b(2)],'r:')
plot([0,T-1],[b(3),b(3)],'g:')
hold off
ylabel('b(t)')
subplot(3,1,3)
plot(time,SS2)
ylabel('s^2(t)')
xlabel('time')







