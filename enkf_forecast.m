function [xf] = enkf_forecast(model, t0, t1, xa)
%ENKF_FORECAST EnKF forecast step.
%   XF = ENKF_FORECAST(MODEL, T0, T1, XA)
%   See Software Design Notes for more details.

% $Id: enkf_forecast.m,v 1.2 2005/07/17 13:38:50 dlivings Exp $
[N,n] = size(xa);


xf = zeros(length(t1), N, n);

for i = 1 : n
     xf(:,:,i) = model(t0, t1, xa(:,i));
end
