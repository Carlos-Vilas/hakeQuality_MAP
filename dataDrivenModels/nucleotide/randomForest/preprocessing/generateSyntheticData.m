function Xsynth = generateSyntheticData(X, nRep)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%       
% Function to generate synthetic data from the experimental data
% Input arguments:
% - X:     Real data matrix (rows = times, columns = repetitions).
% - nRep:  Number of synthetic repetitions to generate per time point.
% Output argument:
% - Xsynth: Generated synthetic data (rows = times, columns = nRep)
%           One synthetic sample per time x repetition.

% Number of time points
[nT, ~] = size(X);

% Output synthetic matrix
Xsynth = NaN(nT, nRep);

% Loop through each time point
for ii = 1:nT

    % Extract repetitions for this time (remove NaNs)
    x = X(ii, :);
    x = x(~isnan(x));

    % If there are no valid observations
    if isempty(x)
        Xsynth(ii, :) = NaN;
        continue;
    end

    % Admissible margins beyond observed values
    delta = 0.6;                       % precentage
    a     = min(0,min(x)*(1.0-delta)); % minimum
    b     = max(x)*(1.0+delta);        % maximum

    % Prior hyperparameters
    mu0    = 0;
    kappa0 = 0.05;
    alpha0 = 2;
    beta0  = 1;

    % Features
    n    = numel(x);
    xbar = mean(x);

    % Posterior hyperparameters
    kappa_n = kappa0 + n;
    mu_n    = (kappa0*mu0 + n*xbar) / kappa_n;
    alpha_n = alpha0 + n/2;
    if n > 1
        s2 = var(x,1);
    else
        s2 = 0;
    end
    beta_n = beta0 + 0.5*n*s2 + (kappa0*n*(xbar - mu0)^2)/(2*kappa_n);

    % Sample sigma^2 and mu
    sigma2 = 1 ./ gamrnd(alpha_n, 1./beta_n);
    mu     = mu_n + sqrt(sigma2/kappa_n)*randn;

    % Standardized truncation limits
    alpha_t = (a - mu)/sqrt(sigma2);
    beta_t  = (b - mu)/sqrt(sigma2);

    % Truncated normal sampling using inverse CDF
    u  = rand(nRep,1);
    Fa = normcdf(alpha_t);
    Fb = normcdf(beta_t);
    z  = norminv(Fa + u.*(Fb - Fa));

    % Final synthetic data
    xsyn = mu + sqrt(sigma2)*z;

    % Add the data to the current time point
    Xsynth(ii, :) = xsyn;
end

% End function
end