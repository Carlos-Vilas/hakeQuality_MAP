function J = costWrapper_06(eZ, otherData)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%    
% Function to call the model simulation program with the data of each
% experiment included in the parameter estimation and compute the Gaussian
% negative log-likelihood
% Input arguments:
% - eZ:        Unknown model parameters in log scale
% - otherData: Structure containing the experimental data, and the 
%              known model parameters
% Output argument:
% - J:  Negative log-likelihood considering all the experiments used in the
%       parameter estimation procedure     

% Number of physical parameters
npp = otherData.npp;

% Convert parameters from log-space to physical space
eP_model = exp(eZ(1:npp));

% One scale per variable (all strictly positive)
s_Sh  = exp(eZ(npp+1));
s_Ph  = exp(eZ(npp+2));
s_TVB = exp(eZ(npp+3));
s_IMP = exp(eZ(npp+4));
s_Ino = exp(eZ(npp+5));
s_Hx  = exp(eZ(npp+6));

% Parameter-domain check
if any(eP_model <= 0.0) || any(~isfinite([s_Sh s_Ph s_TVB s_IMP s_Ino s_Hx]))
    J = 1.0e10;
    return
end

% Unpack data
kP      = otherData.kP;
expData = otherData.expData;
nExp    = length(expData);

% Preallocation
Jexp = zeros(1, nExp);
tSim = cell(1, nExp);
ySim = cell(1, nExp);

% Degrees of freedom for Student's t
nu      = otherData.nu;

% Model output indices
iSh = 1;  iPh = 2;  iTVB = 3;  iIMP = 4;  iIno = 5;  iHx = 6;

% Model simulation and NLL computation
try
    for iExp = 1:nExp
    
        % Simulate experiment iExp
        [tSim{iExp}, ySim{iExp}] = simulationModel06(eP_model, kP,...
                expData{iExp}.tT, expData{iExp}.T, expData{iExp}.sT,...
                iExp, nExp);
        
        % Check whether any of the values in ySim is NaN (failed
        % simulation)
        nanCheck = isnan(ySim{iExp});
        flagNan  = sum(nanCheck(:));

        % Computation of the NLL
        if flagNan ~= 0 
            Jexp(iExp) = 1.0e10;
        else 
            % Interpolate model outputs at sampling times
            ySh_ts  = interp1(tSim{iExp}, ySim{iExp}(:,iSh),  expData{iExp}.tsSh);
            yPh_ts  = interp1(tSim{iExp}, ySim{iExp}(:,iPh),  expData{iExp}.tsPh);
            yTVB_ts = interp1(tSim{iExp}, ySim{iExp}(:,iTVB), expData{iExp}.tsTVB);
            yIMP_ts = interp1(tSim{iExp}, ySim{iExp}(:,iIMP), expData{iExp}.tsIMP);
            yIno_ts = interp1(tSim{iExp}, ySim{iExp}(:,iIno), expData{iExp}.tsIno);
            yHx_ts  = interp1(tSim{iExp}, ySim{iExp}(:,iHx),  expData{iExp}.tsHx);

            % Student-t NLL for all repetitions
            nllSh  = computeNLL_reps(ySh_ts,  expData{iExp}.Sh,  s_Sh,  nu);
            nllPh  = computeNLL_reps(yPh_ts,  expData{iExp}.Ph,  s_Ph,  nu);
            nllTVB = computeNLL_reps(yTVB_ts, expData{iExp}.TVB, s_TVB, nu);
            nllIMP = computeNLL_reps(yIMP_ts, expData{iExp}.IMP, s_IMP, nu);
            nllIno = computeNLL_reps(yIno_ts, expData{iExp}.Ino, s_Ino, nu);
            nllHx  = computeNLL_reps(yHx_ts,  expData{iExp}.Hx,  s_Hx,  nu);

            Jexp(iExp) = nllSh + nllPh + nllTVB + nllIMP + nllIno + nllHx;
        end
    end

catch
    % Penalty if simulation fails
    J = 1.0e10;
    return
end

% Total cost
J = sum(Jexp);

if ~isfinite(J) || ~isreal(J)
    J = 1.0e10;
    return
end

end