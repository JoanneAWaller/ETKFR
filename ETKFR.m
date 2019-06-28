function  [dob,doa,MEAN_XF,MEAN_XA,EstR] = ETKFR(model,t,kl,xf0,y,H,AllR,Ns,regularisation)
tic
% Determine the total number of time steps
kfinal = length(t);
% Determine the number of gridpoints and ensembles
XF = xf0;

Nm = size(xf0,1);
Np = size(y,2);

% Determine the total number of assimilation steps
lfinal = length(kl);

%Determine storage matrices
MEAN_XF = zeros(Nm,lfinal);
MEAN_XA = zeros(Nm,lfinal);
dob = zeros(lfinal,Np);
doa = zeros(lfinal,Np);
EstR = zeros(Np,Np,kfinal);

%% Run model until given time
if lfinal == 0
    % If there are no observations then forecast to the final time.
    k = 2:kfinal;
    [F] = enkf_forecast(model, t(1), t(k), XF);
    b = length(k);
    F = squeeze(F(b,:,:));
    XF = F;
elseif kl(1) > 1
    % Forecast to first observation.
    k = 2:kl(1);
    [F] = enkf_forecast(model, t(1), t(k), XF);
    b = length(k);
    F = squeeze(F(b,:,:));
    XF = F;
end

%% Assimilation steps
waitbar(0,'Please wait. Runing assimilation');
for l = 1 : lfinal
    waitbar(l/lfinal)
    
    % Set R
    if l <=Ns
        R = squeeze(AllR(:,:,l));
        rR = sqrtm(R);
    else
        % Estimate R (If using the ETKFR)
        as = doa(l-Ns:l-1,:);
        bs = dob(l-Ns:l-1,:);
        R = as'*bs/(Ns-1);
        Rsym = 0.5*(R + R');
        
        [R] = regularisation(Rsym); 
        EstR(:,:,l) = R;
        
        rR = sqrtm(R); 
    end
    
    
    
    %Assimilate next observation.
    [XA(:,:),MEAN_XF(:,l),MEAN_XA(:,l)] = etkf_analysis(XF, y(l,:)', H, rR);
        
    % Store the background innovations
    yf = H*XF;
    mean_yf = mean(yf, 2);
    dob(l,:) = y(l,:)' - mean_yf;
     
    % Store the analysis innovations
    ya = H*XA;
    mean_ya = mean(ya, 2);
    doa(l,:) = y(l,:)' -  mean_ya;
        
    
    %% Forecast ensemble members
    if l < lfinal
        % More observations coming up; forecast to next one.
        k = kl(l)+1 : kl(l+1);
        [F] = enkf_forecast(model, t(kl(l)), t(k), XA);
        b = length(k);
        F = squeeze(F(b,:,:));
        XF = F;
    elseif kl(l) < kfinal;
        % Last observation; forecast to end.
        k = kl(l)+1 : kfinal;
        [F] = enkf_forecast(model, t(kl(l)), t(k), XA);
        b = length(k);
        F = squeeze(F(b,:,:));
        XF = F;
    end
end
toc
