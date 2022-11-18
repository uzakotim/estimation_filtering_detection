% Linear
% Exponential
% Time Varying Exponential
% Restricted
% Use of prior information
function theta = forgetting(input,output,A0,B0,n_samples,n_samplesBatch,variance)
    a = A0;
    b = B0;
    na = numel(a)-1;
    nb = numel(b);
    f = @(z, a, b) a(2:end)*z(1:na) + b*z(na + 1:end) + sqrt(variance)*randn();
    
    y_est = zeros(n_samples, 1);
    y_est(4:n_samplesBatch) = output(4:n_samplesBatch);
    y = output;
    u = input;
    theta = zeros(na + nb, n_samples);
    P=sqrt(variance)*eye(7,7);
    for i = na+1:n_samples
        % Regresor construction
        z = [-y(i - 1:-1:i-na); u(i:-1:i-nb + 1)];
        y_est(i) = f(z, a, b);
        
        % Zeta
        zeta = z'*P*z;
        
        % Prediction error
        epsilon = y(i) - z'*theta(:, i-1);
        
        % Kalman gain
        K = P*z/(1+zeta);
        
        %%  Filtering (Data) step:
        %   Mean value update
        theta(:,i) = theta(:,i-1) + K*epsilon;
    
        %   Covariance update
        P = P - P*z*z'*P/(1+zeta);
        
        P_data_step = P;
        
        %%  Propagation (Time) step:
        %   mean value parameter estimation does not change
        %       theta(t+1|t) = theta(t|t)
        %
        %   updating the covariance matrix
        
        %%%  No forgetting
%              P = P;
%         
        %%%  Linear forgetting
    %         V = 0.000001*eye(na + nb);
    %         P = P + V;
    %     
        %%%  Exponential forgetting
%             tef = 1500;
%             phi=1-1/tef;
%             P=1/phi*P;
%         
        %%%  Restricted linear forgetting
            K = P*z;
            V = 0.00000075*eye(7);
            zeta_v = z'*V*z;
            Vo = K*(zeta_v/zeta^2)*K';
            P = P + Vo;
        
        %%%  Restricted exponential forgetting
%             K = P*z;
%             tef = 250;
%             fi=1-1/tef;
%             zeta_v = (1-fi)/fi * zeta;
%             Vo = K*(zeta_v/zeta^2)*K';
%             P = P + Vo;
%         
        P_time_step = P;
    end
    return;
end