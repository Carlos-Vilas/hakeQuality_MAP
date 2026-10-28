% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Main file to estimate the parameters of model 03 that considers:
% - Baranyi and Roberts model for bacteria
% - Luedeking-Piret(both bacteria: Sh, Ph) formation of TVB-N. 
% - Nucleotide reactions: r1, r2, r3, r4, r1bac, r2bac, D. No bacterial
%                         catalysis
% Temperature is expressed in Celsius, time in days
% Reparametrized version: parameters are optimized in log-space (z = log(p))
clear variables
%clearvars -except modCon i
clc

% Load and organize experimental data
cd ../../../experimentalData
iComb   = 2;
purpose = 'Est';
expData = readExperimentalData(iComb,purpose);
cd ../knowledgeDrivenModels/bads/model03

% Parallel pool (process-based)
nExp   = length(expData);
nCores = 4;

pp = gcp('nocreate');
if nExp > 1
    if isempty(pp)
        parpool('local', nCores); 
    end
    p = genpath(pwd);
    pctRunOnAll addpath(p);
end

% Known parameters
T0Sh = -9.0;     % Temperature at which there is no growth [Celsius]. From Safety and quality issues in fish processing book Table 12.4. Pag 209
T0Ph = -9.0;     % Temperature at which there is no growth [Celsius]. From Safety and quality issues in fish processing book Table 12.4. Pag 209
A1   = 1.15e20;  % Pre-exponential coefficient in Arrhenius expression (reaction r1) [umol/(d·g-fish)]
A2   = 7.63e15;  % Pre-exponential coefficient in Arrhenius expression (reaction r2) [umol/(d·g-fish)]
A3   = 2.29e18;  % Pre-exponential coefficient in Arrhenius expression (reaction r3) [umol/(d·g-fish)]
A4   = 4.84e12;  % Pre-exponential coefficient in Arrhenius expression (reaction r4) [umol/(d·g-fish)]
R    = 8.314;    % Universal gas constant [J/(mol·K)]
kP   = [T0Sh, T0Ph, A1, A2, A3, A4, R];

% Bundle data for wrappers
otherData.expData = expData;
otherData.kP      = kP;
otherData.nu      = 4; % Degrees of freedom for Student's t likelihood

% Physical parameter bounds
% Bounds on the initial conditions
for iExp = 1:nExp
    
    Sh0_Ub(1,iExp) = 4.0;   % Shewanella [log10(CFU/g-fish)]
    Sh0_U(1,iExp)  = 3.0;   % Shewanella [log10(CFU/g-fish)]
    Sh0_L(1,iExp)  = 1.0;   % Shewanella [log10(CFU/g-fish)]
    Sh0_Lb(1,iExp) = 0.0;   % Shewanella [log10(CFU/g-fish)]
    Sh0(1,iExp)    = 2.0;   % Shewanella [log10(CFU/g-fish)]

    Ph0_Ub(1,iExp) = 6.0;   % Photobacterium [log10(CFU/g-fish)]
    Ph0_U(1,iExp)  = 5.9;   % Photobacterium [log10(CFU/g-fish)]
    Ph0_L(1,iExp)  = 2.5;   % Photobacterium [log10(CFU/g-fish)]
    Ph0_Lb(1,iExp) = 2.0;   % Photobacterium [log10(CFU/g-fish)]
    Ph0(1,iExp)    = 3.0;   % Photobacterium [log10(CFU/g-fish)]

    TVB0_Ub(1,iExp) = 20.0/100.0;  % Total volatile base nitrogen [g/g-fish]             
    TVB0_U(1,iExp)  = 18.0/100.0;  % Total volatile base nitrogen [g/g-fish]             
    TVB0_L(1,iExp)  = 7.0/100.0;   % Total volatile base nitrogen [g/g-fish]             
    TVB0_Lb(1,iExp) = 4.0/100.0;   % Total volatile base nitrogen [g/g-fish]             
    TVB0(1,iExp)    = 10.0/100.0;  % Total volatile base nitrogen [g/g-fish]             

    IMP0_Ub(1,iExp) = 3.0;  % Inosine monophosphate [umol/g-fish]
    IMP0_U(1,iExp)  = 2.5;  % Inosine monophosphate [umol/g-fish]
    IMP0_L(1,iExp)  = 1.0;  % Inosine monophosphate [umol/g-fish]
    IMP0_Lb(1,iExp) = 0.3;  % Inosine monophosphate [umol/g-fish]
    IMP0(1,iExp)    = 2.0;  % Inosine monophosphate [umol/g-fish]

    Ino0_Ub(1,iExp) = 0.9;  % Inosine [umol/g-fish]
    Ino0_U(1,iExp)  = 0.6;  % Inosine [umol/g-fish]
    Ino0_L(1,iExp)  = 0.05; % Inosine [umol/g-fish]
    Ino0_Lb(1,iExp) = 0.0;  % Inosine [umol/g-fish]
    Ino0(1,iExp)    = 0.1;  % Inosine [umol/g-fish]

    Hx0_Ub(1,iExp) = 0.5;   % Hypoxanthine [umol/g-fish]
    Hx0_U(1,iExp)  = 0.3;   % Hypoxanthine [umol/g-fish]
    Hx0_L(1,iExp)  = 0.05;  % Hypoxanthine [umol/g-fish]
    Hx0_Lb(1,iExp) = 0.0;   % Hypoxanthine [umol/g-fish]
    Hx0(1,iExp)    = 0.1;   % Hypoxanthine [umol/g-fish]

end
% Bounds on the global parameters
%      bSh,   Shmax, bPh,   Phmax, YSh,   mSh,   YPh,   mPh,   Ea1,   Ea2,   Ea3,   Ea4,   D 
ub  = [1.2e0, 8.5,   1.2e1, 8.5,   1.0e3, 1.5e3, 1.0e3, 1.5e3, 4.0e5, 4.0e5, 4.0e5, 4.0e5, 1.0,...
       Sh0_Ub, Ph0_Ub, TVB0_Ub, IMP0_Ub, Ino0_Ub, Hx0_Ub];

%      bSh,   Shmax, bPh,   Phmax, YSh,   mSh,    YPh,   mPh,   Ea1,   Ea2,   Ea3,   Ea4,   D 
pub = [1.0e0, 8.0,   1.0e1, 8.0,   0.9e3, 1.4e3, 0.9e3, 1.4e3, 3.0e5, 3.0e5, 3.9e5, 3.9e5, 0.9,...
       Sh0_U, Ph0_U, TVB0_U, IMP0_U, Ino0_U, Hx0_U];

%      bSh,     Shmax, bPh,    Phmax,  YSh,     mSh,     YPh,     mPh,     Ea1,   Ea2,   Ea3,   Ea4,   D 
plb = [2.0e-12, 4.0,   2.0e-8, 5.0,    1.0e-11, 1.0e-11, 1.0e-11, 1.0e-11, 4.0e4, 4.0e4, 4.0e4, 4.0e4, 1.0e-13,...
       Sh0_L, Ph0_L, TVB0_L, IMP0_L, Ino0_L, Hx0_L];

%      bSh,     Shmax, bPh,    Phmax,  YSh,     mSh,     YPh,     mPh,     Ea1,   Ea2,   Ea3,   Ea4,   D 
lb  = [1.0e-12, 3.9,   1.0e-8, 4.9,    1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 3.5e4, 3.5e4, 3.5e4, 3.5e4, 1.0e-14,...
       Sh0_Lb, Ph0_Lb, TVB0_Lb, IMP0_Lb, Ino0_Lb, Hx0_Lb];

% Number of physical parameters
otherData.npp = length(ub); 

% Extend the parameters to be estimated to include sigmas
%           Shs     Phs     TVBs    IMPs    Inos    Hpxs
ub  = [ub,  10.0,   10.0,   10.0,   10.0,   10.0,   10.0];
pub = [pub, 1.0,    1.0,    1.0,    1.0,    1.0,    1.0];
plb = [plb, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4];
lb  = [lb,  1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6];

% Number of runs
nRuns = 100;  

% Ensure parameters ar possitve so they can be transformed to the log-space
eps_pos = 1e-14;
lb_pos  = max(lb,  eps_pos);
ub_pos  = max(ub,  eps_pos);
plb_pos = max(plb, eps_pos);
pub_pos = max(pub, eps_pos);
% Convert all to log-space
z_lb    = log(lb_pos);
z_ub    = log(ub_pos);
z_plb   = log(plb_pos);
z_pub   = log(pub_pos);

% Number of parameters (same dimensionality as your parameter vector)
nParams = numel(z_plb);

% Preallocate and sample LHS in [0,1]
z0 = zeros(nRuns, nParams);
lhsUnit = lhsdesign(nRuns, nParams);

% Scale each dimension to [z_plb(ii), z_pub(ii)]
for ii = 1:nParams
    % Linear mapping from [0,1] to [z_plb(ii), z_pub(ii)]
    z0(:, ii) = z_plb(ii) + lhsUnit(:, ii) .* (z_pub(ii) - z_plb(ii));
end

% BADS options
options = bads('defaults');
options.TolMesh             = 1e-6;
options.MaxFunEvals         = 1.0e5;
options.UncertaintyHandling = false;

% Call the bads toolbox
p_est = zeros(nRuns, nParams);
J_est = zeros(nRuns, 1);
parfor ii = 1:nRuns
    % Call bads
    [z_est_i, J_est_i] = bads(@(z) costWrapper_03(z, otherData), ...
                 z0(ii,:), z_lb, z_ub, z_plb, z_pub, [], options);

    % Store results in p space
    p_est(ii,:) = exp(z_est_i); 
    J_est(ii)   = J_est_i; 
end

% Remove previously created paths
pctRunOnAll rmpath(p);

% Minimum value of J_est
[Jbest, idx] = min(J_est);
Pbest        = p_est(idx,:);

% Save the results
save PEresults_model03 Pbest Jbest J_est p_est

