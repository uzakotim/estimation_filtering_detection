% likelihood of sigma - normal distribution, known mu
s2_statistics = 1:20;
xx = 0.2:0.1:100;
yy = xx;
N = 5;
n=length(xx);
figure(1)
hold off
for s = s2_statistics,
    for k=1:n,
        yy(k)=(1/xx(k))^N * exp(-N*s^2/(2*xx(k)^2));
    end
    yy = yy/sum(yy);
    plot(xx,yy)
    hold on
    grid on
end
plot(xx,.12./xx,'r','linewidth',2)
xlim([0, 50])
ylim([0,0.15])
title('Likelihood l(\sigma|s^2 = 1:20) vs. 1/\sigma')
xlabel('variance \sigma')
%
figure(2)
hold off
for s = s2_statistics,
    for k=1:n,
        yy(k)=exp( N*(log(s)-log(xx(k))) - N/2*exp(2*(log(s)-log(xx(k)))));
    end
    plot(log(xx),yy)
    hold on
    grid on
end
%plot(xx,.12./xx,'r','linewidth',2)
xlim([-1, 5])
% ylim([0,0.09])
title('Likelihood l(log(\sigma)|s^2 = 1:20)')
xlabel('log(\sigma)')
