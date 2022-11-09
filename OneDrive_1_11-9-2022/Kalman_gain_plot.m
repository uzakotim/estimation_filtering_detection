% MS/LMS estimate for normal distribution
%
P = [2 1;1 0.51];
sqrtP = chol(P)';
%
phi=0:0.02:6.4;
y1 = cos(phi);
y2 = sin(phi);
y = [y1;y2];
%
% plot(y1,y2);
%
% contours for alpha = 1, 2, 4
mu = [1;1];
x_1 = mu+sqrtP*y*sqrt(1);
x_2 = mu+sqrtP*y*sqrt(2);
x_4 = mu+sqrtP*y*sqrt(4);
%
% figure(1)
plot(x_1(1,:),x_1(2,:))
axis('equal'); grid on
hold on
plot(x_2(1,:),x_2(2,:))
plot(x_4(1,:),x_4(2,:))
title('PDF contours / LMS gain'); xlabel('x'); ylabel('y');
%
% semiaxes
%
[V,D] = eig(P);
v1 = V(:,1)*sqrt(D(1,1));
v2 = V(:,2)*sqrt(D(2,2));
plot([mu(1), mu(1)+v1(1)],[mu(2), mu(2)+v1(2)],'r','linewidth',2)
plot([mu(1), mu(1)-v1(1)],[mu(2), mu(2)-v1(2)],'r','linewidth',2)
plot([mu(1), mu(1)+v2(1)],[mu(2), mu(2)+v2(2)],'r','linewidth',2)
plot([mu(1), mu(1)-v2(1)],[mu(2), mu(2)-v2(2)],'r','linewidth',2)
%
% Kalman gain
%
K = P(1,2)/P(2,2);
dY=-3:3;
dX=K*dY;
plot(mu(1)+dX,mu(2)+dY,'b','linewidth',2)
% hold off

