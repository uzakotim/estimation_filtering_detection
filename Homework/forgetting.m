% Linear
% Exponential
% Time Varying Exponential
% Restricted
% Use of prior information
function [] = forgetting(input,output,initCov,ForgettingFactor,A0,B0,a1_,a2_,a3_,b1_,b2_,b3_,ttl)
    na = 3;
    nb = 4;
    nk = 0;
    method = 'ForgettingFactor';
    obj = recursiveARX([na nb nk],A0,B0,'EstimationMethod',method);
    obj.InitialParameterCovariance = initCov;
    obj.ForgettingFactor = ForgettingFactor;
    
    A = zeros(numel(input),4);
    B = zeros(numel(input),4);
    EstimatedOutput = zeros(numel(input),1);
    for i = 1:numel(input)
        [A(i,:),B(i,:),EstimatedOutput(i,:)] = step(obj,output(i),input(i));
    end
    
    subplot(2,1,1)
    plot(1:1:length(input),A);
    hold on
    plot([0,1000,1000,2000,2000,3000],[a1_(1,1),a1_(1,1),a2_(1,1),a2_(1,1),a3_(1,1),a3_(1,1)],'b:')
    plot([0,1000,1000,2000,2000,3000],[a1_(1,2),a1_(1,2),a2_(1,2),a2_(1,2),a3_(1,2),a3_(1,2)],'b:')
    plot([0,1000,1000,2000,2000,3000],[a1_(1,3),a1_(1,3),a2_(1,3),a2_(1,3),a3_(1,3),a3_(1,3)],'b:')
    plot([0,1000,1000,2000,2000,3000],[a1_(1,4),a1_(1,4),a2_(1,4),a2_(1,4),a3_(1,4),a3_(1,4)],'b:')
    hold off
    ylabel('a(t)')
    title(ttl);
    
    subplot(2,1,2)
    plot(1:1:length(input),B)
    hold on
    plot([0,1000,1000,2000,2000,3000],[b1_(1,1),b1_(1,1),b2_(1,1),b2_(1,1),b3_(1,1),b3_(1,1)],'b:')
    plot([0,1000,1000,2000,2000,3000],[b1_(1,2),b1_(1,2),b2_(1,2),b2_(1,2),b3_(1,2),b3_(1,2)],'b:')
    plot([0,1000,1000,2000,2000,3000],[b1_(1,3),b1_(1,3),b2_(1,3),b2_(1,3),b3_(1,3),b3_(1,3)],'b:')
    plot([0,1000,1000,2000,2000,3000],[b1_(1,4),b1_(1,4),b2_(1,4),b2_(1,4),b3_(1,4),b3_(1,4)],'b:')
    hold off 
    ylabel('b(t)')
end