function [stheta_out, ssigma_out, k, eps, dy] = ldfil(stheta, ssigma, data, phi)
% =========================================================
% LDFIL - ARX identification with LD factorized covariance
% =========================================================
% stheta = [theta,d,L^T]
% ssigma = [nu,nus2]
% data = [z^T,y]
% phi = forgetting coefficient
% =========================================================
% unpack input ============================================
n = length(data) - 1;
z = data(1:n)';
y = data(n+1);
theta = stheta(:,1);
dth = stheta(:,2);
ltht = stheta(:,3:n+2);
% check forgetting factor =================================
if nargin < 4
    phi = 1;
else
    phi = min(phi,1);
    phi = max(phi,.01);
end
% update covariance =======================================
m = [1, zeros(1,n); ltht*z, ltht ];
d = [1; dth ];
% run dydr for i=1 ========================
for j = n+1:-1:2
    mji = m(j,1);
    di = d(1)+mji*mji*d(j);
    mu = mji*d(j)/di;
    d(j) = d(j)*d(1)/di;
    d(1) = di;
    m(j,:) = m(j,:)-mji*m(1,:);
    m(1,:) = m(1,:)+mu*m(j,:);
end
% decompose ===============================
dy = d(1);
dth = d(2:n+1)/phi;
ltht = m(2:n+1,2:n+1);
k = m(1,2:n+1)';
% update theta ============================
eps = y - z'*theta;
theta = theta + k*eps;
stheta_out = [theta,dth,ltht];
% update sigma ============================
if ~isempty(ssigma)
    nu = ssigma(1);
    nus2 = ssigma(2);
    nu = phi*(nu+1);
    nus2 = phi*(nus2+eps*eps/dy);
    ssigma_out = [nu,nus2];
else
    ssigma_out = [];
end
end
