function logpost = costWrapper_07_logJoint(eZ, otherData, prior_z)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Function to compute the Log-posterior for VBMC: Student-t log-likelihood
% + LogNormal prior. 
% Input arguments:
% - eZ:        Unknown model parameters in log scale
% - otherData: Structure containing the experimental data, and the 
%              known model parameters
% - prior_z:   Prior selected by the user (the value obtained in BADS is
%              used)
% Output argument:
% - logpost:  Log-posterior = Student-t log-likelihood + LogNormal prior.    
% Parameters are given in log-space as eZ. We convert to physical space eP = exp(eZ).

% Convert parameters to physical space
eP = exp(eZ);

% Penalization for invalid points
logpost_bad = -1.0e5 - sum((eZ - prior_z).^2.0);

% Parameter-domain check (positivity)
if any(eP <= 0.0) || any(~isfinite(eP))
    logpost = logpost_bad;
    return
end

% Unpack data
kP      = otherData.kP;
expData = otherData.expData;
nExp    = length(expData);

% Student-t degrees of freedom
nu = otherData.nu;

% Number of physical parameters
npp = otherData.npp;

% Extract physical parameters
eP_model = eP(1:npp);

% Extract scales (one per variable)
s_Sh  = eP(npp+1);
s_Ph  = eP(npp+2);
s_TVB = eP(npp+3);
s_IMP = eP(npp+4);
s_Ino = eP(npp+5);
s_Hx  = eP(npp+6);

% Invalid sigma
if any([s_Sh s_Ph s_TVB s_IMP s_Ino s_Hx] <= 0.0) || ...
   any(~isfinite([s_Sh s_Ph s_TVB s_IMP s_Ino s_Hx]))
    logpost = logpost_bad;
    return
end

% Preallocation
loglik_exp = zeros(1, nExp);

% Model output indices
iSh = 1; iPh = 2; iTVB = 3; iIMP = 4; iIno = 5; iHx = 6;

try
    for iExp = 1:nExp

        % Simulate experiment
        [tSim{iExp}, ySim{iExp}] = simulationModel07(eP_model, kP, ...
            expData{iExp}.tT, expData{iExp}.T, expData{iExp}.sT, ...
            iExp, nExp );

        % Interpolate model at sampling times
        ySh_ts  = interp1(tSim{iExp}, ySim{iExp}(:,iSh),  expData{iExp}.tsSh);
        yPh_ts  = interp1(tSim{iExp}, ySim{iExp}(:,iPh),  expData{iExp}.tsPh);
        yTVB_ts = interp1(tSim{iExp}, ySim{iExp}(:,iTVB), expData{iExp}.tsTVB);
        yIMP_ts = interp1(tSim{iExp}, ySim{iExp}(:,iIMP), expData{iExp}.tsIMP);
        yIno_ts = interp1(tSim{iExp}, ySim{iExp}(:,iIno), expData{iExp}.tsIno);
        yHx_ts  = interp1(tSim{iExp}, ySim{iExp}(:,iHx),  expData{iExp}.tsHx);

        % Student‑t log-likelihood = − NLL
        llkSh  = -computeNLL_reps(ySh_ts,  expData{iExp}.Sh,  s_Sh,  nu);
        llkPh  = -computeNLL_reps(yPh_ts,  expData{iExp}.Ph,  s_Ph,  nu);
        llkTVB = -computeNLL_reps(yTVB_ts,expData{iExp}.TVB, s_TVB, nu);
        llkIMP = -computeNLL_reps(yIMP_ts,expData{iExp}.IMP, s_IMP, nu);
        llkIno = -computeNLL_reps(yIno_ts,expData{iExp}.Ino, s_Ino, nu);
        llkHx  = -computeNLL_reps(yHx_ts, expData{iExp}.Hx,  s_Hx,  nu);

        loglik_exp(iExp) = llkSh + llkPh + llkTVB + llkIMP + llkIno + llkHx;
    end

catch
    logpost = logpost_bad;
    return
end

% Total log-likelihood
loglik = sum(loglik_exp);

if ~isfinite(loglik) || ~isreal(loglik)
    logpost = logpost_bad;
    return
end

% Log-normal prior in log-space
cv = 0.5;
logprior_terms = -0.5*((eZ - prior_z)./cv).^2 - log(cv) - 0.5*log(2*pi);
logprior = sum(logprior_terms);

if ~isfinite(logprior)
    logpost = logpost_bad;
    return
end

% Final log-posterior
logpost = loglik + logprior;

if ~isfinite(logpost) || ~isreal(logpost)
    logpost = logpost_bad;
end

% Clip extremely large values for VBMC stability
logpost = max(min(logpost, 1e6), -1.0e6);

end
