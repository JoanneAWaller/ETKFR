function [xa,mean_xf,mean_xa] = etkf_analysis(xf, y, H, rR)
%ETKF analysis step

n = size(xf,2);

% Determine the ensemble mean
mean_xf = mean(xf, 2);
% Determine the ensembel perterbation matrix
Xpf = (xf - mean_xf*ones(1, n)) / sqrt(n-1);
[b,c] = size(Xpf);

% Determine the B matrix
B = Xpf*Xpf';

% Detemine the matrix holding the mesuremetns of the perturbations
yf = H*xf;
mean_yf = mean(yf, 2);
Ysf = H*Xpf;
Ysf = rR \ Ysf;

[U, S, V] = svd(Ysf');
XpfU = Xpf*U;

A = eye(c) + S*S';

Xpa = XpfU / sqrtm(A) * U;

g = 1;

d = (y - mean_yf);
d = rR \ d;
d = V'*d;
d = (S'*S + eye(length(y))) \ d;
d = S*d;
d = XpfU*d;
mean_xa = mean_xf + d;

xa = mean_xa*ones(1, n) + g*(sqrt(n-1)*Xpa - (1/n)*sqrt(n-1)*Xpa*ones(n,1)*ones(1,n));