function Jnll = computeNLL_reps(yts, yMeas, s, nu)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Compute Student-t NLL using all repetitions, ignoring NaNs
% Input arguments:
% - yts:   Column vector of theoretical values (nTimes x 1)
%          conditions (yTVB0 in each of the experiements used for estimation)
% - yMeas: Matrix of measurements (nTimes x nReps)
% - nu:    Degrees of freedom of the Student-t distribution
% - s:     Scale parameter of the Student-t distribution
% Output argument:
% - Jnll:  Negative log-likelihood 

% Replicate yts to match yMeas size
Y = repmat(yts, 1, size(yMeas,2));   % nTimes x nReps

% Constant term in Jnll
const_t = 0.5*log(pi*nu) + gammaln(nu/2) - gammaln((nu+1)/2);

% Compute residuals
R = Y - yMeas;       

% Ignore NaNs
mask = ~isnan(R);

% Vectorized NLL contribution
term = log(s) + const_t + 1.0/2.0*(nu+1.0)*log(1.0 + ((R./s).^2.0)/nu);

% Apply mask and sum
Jnll = sum(term(mask));

end