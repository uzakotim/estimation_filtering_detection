close all
%% Filtering Circuit
% The physical parameters for our example are:
% (L1)   inductance from first loop     2 H
% (C1)   capacitance from first loop    0.05 F
% (L2)   inductance from second loop    1 H
% (C2)   capacitance from second loop   0.047 F
% (R)    resistance                     10 Ohm

L1 = 2;
C1 = 0.05;
L2 = 1;
C2 = 0.047;
R = 10;

% e(t) modeled as white random process e(t) = w(t)

% x(t) = [q1(t) q2(t) q1(t)' q2(t)']
% q1(t)'' = 1/L1*(-R(q1(t)'-q2(t)')-1/C1*q1(t))
% q2(t)'' = 1/L2*(-R(q2(t)'-q1(t)')-1/C2*q2(t))
% xdot(t) = Ac*x(t) + G*w(t), w(t) ∼ N(0,Qc) 
% y(t) = Cc*x(t) + e(t), e(t) ~ N(0,Rc) 

Ac = [0 0 1 0; 
    0 0 0 1; 
    -1/(L1*C1) 0 -R/L1 R/L1; 
    0 -1/(L2*C2) R/L2 -R/L2];

Gc = [0 0 1 0]';

C = [1 0 0 0;
    0 1 0 0];

Qc = [10];
Qp = zeros(4,4);
Qp(3,3) = Qc;

Ts = 1;
Tc = 0.2;

Rd = eye(2)*0.25;

%% Discretization

sysC = ss(Ac, Gc, C, []);
sysD = c2d(sysC, Tc);

figure();
step(sysC);
hold on
step(sysD);
hold off
grid on
title('Step response');
legend('Continuous', 'Discretized');
xlim([0, 15]);

Ad = sysD.A;
Gd = sysD.B;
Cd = sysD.C;

Q = @(t) expm(Ac*t)*Qp*expm(Ac'*t);
Qd = integral(Q, 0, Tc, 'ArrayValued', 1);



%% Data generation
T = Tc;
sysD = ss(Ad, eye(4), eye(4), [], Tc);
data_len = 2000;
t = 0:T:(data_len*T)-T;
% w = (sqrt(Qd)*randn(length(Qd), data_len)).';
% n = (randn(length(Rd), data_len).'*sqrt(Rd)).';

w = mvnrnd(zeros(4,1), Qd, data_len).';
n = mvnrnd(zeros(2,1), Rd, data_len).';

Xt = lsim(sysD, w, t).';
Yt = C*Xt + n;

x0 = zeros(4, 1);
P0 = eye(4)*100;

%% Kalman algorithm
xhat = zeros(4, data_len);
xhat(:, 1) = x0;
P = cell(1, data_len);
P{1} = P0;
% innovation covariance
Sk = cell(1, floor(data_len/5));
% innovation
v = zeros(2, floor(data_len/5));
% normalised innovation
q = zeros(1, floor(data_len/5));
% error
err = zeros(2, floor(data_len/5));
for i = 1:data_len
    curx = xhat(:, i);
    curP = P{i};
    % measures
%     if(i ~= 0)
%         v(:, i) = Yt(:, i) - C*curx;
%         Sk{i} = Rd + C*P{i}*C';
%         q(i) = v(:, i)'*inv(Sk{i})*v(:, i);
%     end
    if(mod(i, 5) == 1)
        % next output
        cury = Yt(:, i);
        
        idx = floor(i/5)+1;
        v(:, idx) = cury - C*curx;
        Sk{idx} = Rd + C*P{i}*C';
        q(idx) = v(:, idx)'*inv(Sk{idx})*v(:, idx);
        
        % data-update step
        L = curP*C'*inv(C*curP*C' + Rd);
        curx = curx + L*(cury - C*curx);
        curP = curP - L*C*curP;
        xhat(:, i) = curx;
        err(:, idx) = cury - C*curx;
    end
    % time-update step
    if(i ~= data_len)
        xhat(:, i+1) = Ad*curx;
        P{i+1} = Ad*curP*Ad' + Qd;
    end
end
%% Plotting
plot_est(Yt, xhat, data_len);


%% Results
disp('----------------------Ideal Example----------------------')
[valid, p] = test1(v, Sk);
disp('========Test1========')
disp('Valid')
disp(valid)
disp('Percenatge')
disp(p)

valid = test2(q, 2, length(q));
disp('========Test2========')
disp('Valid')
disp(valid)

disp('========Test3========')
[valid, p] = test3(v);
disp('Valid')
disp(valid)
disp('Percentage')
disp(p)

disp('========RMSE========')
RMSE_err_inn = RMSE(v);
RMSE_err = RMSE(err);
disp('Error from innovation')
disp(RMSE_err_inn)
disp('Error from filtered value')
disp(RMSE_err)


%% Model Mismatch
R = 15;
L1 = 5;
C1 = 0.25;
Ac = [0 0 1 0; 
    0 0 0 1; 
    -1/(L1*C1) 0 -R/L1 R/L1; 
    0 -1/(L2*C2) R/L2 -R/L2];

sysC = ss(Ac, Gc, C, []);
sysD = c2d(sysC, Tc);

Amm = sysD.A;

%% Kalman Algorithm - Model mismatch
xhat = zeros(4, data_len);
xhat(:, 1) = x0;
P = cell(1, data_len);
P{1} = P0;
% innovation covariance
Sk = cell(1, floor(data_len/5));
% innovation
v = zeros(2, floor(data_len/5));
% normalised innovation
q = zeros(1, floor(data_len/5));
% error
err = zeros(2, floor(data_len/5));
for i = 1:data_len
    curx = xhat(:, i);
    curP = P{i};
    
    if(mod(i, 5) == 1)
        % next output
        cury = Yt(:, i);
        
        % measures
        idx = floor(i/5)+1;
        v(:, idx) = cury - C*curx;
        Sk{idx} = Rd + C*P{i}*C';
        q(idx) = v(:, idx)'*inv(Sk{idx})*v(:, idx);
        
        % data-update step
        L = curP*C'*inv(C*curP*C' + Rd);
        curx = curx + L*(cury - C*curx);
        curP = curP - L*C*curP;
        xhat(:, i) = curx;
        err(:, idx) = cury - C*curx;
    end
    % time-update step
    if(i ~= data_len)
        xhat(:, i+1) = Amm*curx;
        P{i+1} = Amm*curP*Amm' + Qd;
    end
end
%% Plotting
plot_est(Yt, xhat, data_len);
%% Results
disp('----------------------Model mismatch----------------------')
[valid, p] = test1(v, Sk);
disp('========Test1========')
disp('Valid')
disp(valid)
disp('Percenatge')
disp(p)

valid = test2(q, 2, length(q));
disp('========Test2========')
disp('Valid')
disp(valid)

disp('========Test3========')
[valid, p] = test3(v);
disp('Valid')
disp(valid)
disp('Percentage')
disp(p)

disp('========RMSE========')
RMSE_err_inn = RMSE(v);
RMSE_err = RMSE(err);
disp('Error from innovation')
disp(RMSE_err_inn)
disp('Error from filtered value')
disp(RMSE_err)

%% Noise Mismatch
Rnm = [0.5 0.1; 0.1 0.75];
Qnm = Qd;
Qnm(2, 2) = 0.5;

%% Kalman Algorithm - Noise mismatch
xhat = zeros(4, data_len);
xhat(:, 1) = x0;
P = cell(1, data_len);
P{1} = P0;
% innovation covariance
Sk = cell(1, floor(data_len/5));
% innovation
v = zeros(2, floor(data_len/5));
% normalised innovation
q = zeros(1, floor(data_len/5));
% error
err = zeros(2, floor(data_len/5));
for i = 1:data_len
    curx = xhat(:, i);
    curP = P{i};
    
    if(mod(i, 5) == 1)
        % next output
        cury = Yt(:, i);
        
        % measures
        idx = floor(i/5)+1;
        v(:, idx) = cury - C*curx;
        Sk{idx} = Rnm + C*P{i}*C';
        q(idx) = v(:, idx)'*inv(Sk{idx})*v(:, idx);
        
        % data-update step
        L = curP*C'*inv(C*curP*C' + Rnm);
        curx = curx + L*(cury - C*curx);
        curP = curP - L*C*curP;
        xhat(:, i) = curx;
        err(:, idx) = cury - C*curx;
    end
    % time-update step
    if(i ~= data_len)
        xhat(:, i+1) = Ad*curx;
        P{i+1} = Ad*curP*Ad' + Qnm;
    end
end
%% Plotting
plot_est(Yt, xhat, data_len);
%% Results
disp('----------------------Noise mismatch----------------------')
[valid, p] = test1(v, Sk);
disp('========Test1========')
disp('Valid')
disp(valid)
disp('Percenatge')
disp(p)

valid = test2(q, 2, length(q));
disp('========Test2========')
disp('Valid')
disp(valid)

disp('========Test3========')
[valid, p] = test3(v);
disp('Valid')
disp(valid)
disp('Percentage')
disp(p)

disp('========RMSE========')
RMSE_err_inn = RMSE(v);
RMSE_err = RMSE(err);
disp('Error from innovation')
disp(RMSE_err_inn)
disp('Error from filtered value')
disp(RMSE_err)




function [valid, percentage] = test1(v, S)
    N = length(v);
    inside = 0;
    valid = 0;
    in = NaN(2, N);
    out = NaN(2, N);
    for i = 1:N
        [V, D] = eig(S{i});
        [d, ind] = sort(diag(D));
        D = D(ind, ind);
        V = V(:, ind);
        A = sqrt(5.9915*D(2, 2))*V(:, 2);
        B = sqrt(5.9915*D(1, 1))*V(:, 1);
        a = norm(A);
        b = norm(B);
        f = sqrt(a^2 - b^2);
        F1 = f*V(:, 2);
        F2 = -f*V(:, 2);
        
        dist = norm(F1 - v(:, i)) + norm(F2 - v(:, i));
        if(dist <= 2*a)
            inside = inside + 1;
            in(:, i) = v(:, i);
        else
            out(:, i) = v(:, i);
        end
    end
    
    percentage = inside/N;
    if(percentage >= 0.95)
        valid = 1;
    end
    
    % ellipse
    [U_calc,D_calc,~] = svd(S{i});
    Ncircle = 100; % number of elements in graph
    theta = [0:1/Ncircle:2*pi+1/Ncircle];
    
    unit_circle(1,:) = cos(theta); % x
    unit_circle(2,:) = sin(theta); % y 
    ell = U_calc*sqrt(5.991*D_calc)*unit_circle;
    
    figure();
    scatter(in(1, :), in(2, :), 'b.');
    hold on
    scatter(out(1, :), out(2, :), 'r.');
    plot(ell(1, :), ell(2, :), 'g', 'LineWidth', 1.5);
    hold off
    title('Innovation boundedness test')
    legend('In 95 % confidence', 'Outliers', 'Covariance ellipsoid')
    grid on
end

function valid = test2(q, m, N)
    qm = mean(q);
    r1 = chi2inv(0.025 ,N*m);
    r2 = chi2inv(0.975, N*m);
    if(qm*N >= r1 && qm*N <= r2)
        valid = 1;
    else
        valid = 0;
    end
end

function [valid, p] = test3(v)
    N = length(v);
    r = zeros(1, length(v));
    valid = 0;
    for tau = 1:N
        for k = 1:(N - tau - 1)
            r(tau) = r(tau) + v(k)'*v(k + tau);
        end
        r(tau) = r(tau)/N;
    end
    bound = 2*(1/sqrt(N));
    in = zeros(1, N);
    in(abs(r) <= bound) = 1;
    p = sum(in)/N;
    if(p >= 0.95)
        valid = 1;
    end
    figure();
    plot(1:N, r, 'LineWidth', 1.5)
    grid on
    title('Autocorrelation')
    hold on
    plot([1,N], [bound, bound], 'r--', 'LineWidth', 1.5)
    plot([1,N], -[bound, bound], 'r--', 'LineWidth', 1.5)
    hold off
    legend('Autocorrelation', '95 % confidence bounds');
end
function RMSE_err = RMSE(err)
    SE = err.^2;
    MSE = mean(SE, 2);
    RMSE_err = sqrt(MSE);
end

function plot_est(y, xhat, data_len)
    ts = 0:1:floor(data_len/5)-1;
    tc = 0:0.2:(data_len*0.2)-0.2;
    y = y(:, 1:5:data_len);
    figure();
    subplot(211)
    plot(ts, y(1, :), 'LineWidth', 1.5)
    hold on
    plot(tc, xhat(1, :) , 'LineWidth', 1.5)
    hold off
    grid on
    title('Charge estimation')
    xlim([150 200]);
    legend('Measured', 'Estimated');
    
    subplot(212)
    plot(ts, y(2, :), 'LineWidth', 1.5)
    hold on
    plot(tc, xhat(2, :), 'LineWidth', 1.5)
    hold off
    grid on
    title('Charge rate estimation')
    xlim([150 200]);
    legend('Measured', 'Estimated');
end

