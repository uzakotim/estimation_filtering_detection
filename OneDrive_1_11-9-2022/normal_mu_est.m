% likelihood of mu - normal distribution, known sigma
y_bar = -5:5;
N = 2;
%
xx = -10:0.1:10;
yy = xx;
n=length(xx);
figure(1)
hold off
for sm = y_bar,
    for k=1:n,
        yy(k)=exp(-(xx(k)-sm)^2/(2/N));
    end
    yy = yy/sqrt(2*pi)/N;
    plot(xx,yy)
    hold on
    grid on
end
title('Likelihood l(\mu|sample mean = -5:5) ')
xlabel('mean \mu')