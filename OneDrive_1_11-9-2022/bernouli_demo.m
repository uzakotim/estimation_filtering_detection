%
% Bayesian estimation of Bernouli distribution parameter \theta
%
Theta = 0.5;
N = 2000;    % number of flops
n = 201;     % number of points
%
% prior distribution
p_theta_prior = ones(1,n);
theta_axis = (0:n-1)/(n-1);
%plot(theta_axis, p_theta_prior)
%axis([-0.1, 1.1, 0, 1.1])
%grid on
%
% flip the coin
%
sample_sum = 0;
%
for f = 1:N,
    if rand()<Theta,
        x=1; 
        % likelihood function
        L = theta_axis;
    else
        x=0;
        L = 1 - theta_axis;
    end;
    sample_sum = sample_sum + x;
    %
    p_theta_posterior = L .* p_theta_prior;
    p_theta_posterior = p_theta_posterior / sum(p_theta_posterior)*(n-1);
    p_max = max([p_theta_posterior, p_theta_prior]);
    p_scale = ceil(100*1.1*p_max)/100;
    %
    plot(theta_axis, p_theta_prior,'b','linewidth',3);
    hold on    
    plot(theta_axis, p_theta_posterior,'r','linewidth',1);
    plot(theta_axis, L,'g','linewidth',2);
    plot(sample_sum/f, 0, 'r*','linewidth',4)
    plot(x, 0, 'go','linewidth',4)
    legend('Prior p(\theta|Y_1^{k-1})','Posterior p(\theta|Y_1^k)','Likelihood L(\theta|y_k)','\theta_{MAP} estimate','Current flip y_k')
    axis([-0.05, 1.05, 0.00, p_scale])
    grid on
    xlabel('\theta');
    ylabel('p(\theta)');
    title(['p(\theta|Y_1^k) after k = ',int2str(f),' coin flips'])
    hold off
    grid on
    if f<10
     input('Hit any key to continue ...');
     clc
    else
     pause(0.001)
    end
    % change of coin ##################################
    if f == 200, Theta = 1 - Theta; end
    %
    p_theta_prior = p_theta_posterior;
end


