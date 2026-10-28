% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Main file to estimate parameter confidence intervals for model 07 that considers:
% - Baranyi and Roberts model for bacteria
% - Power-Law (only Ph) formation of TVB-N. 
% - Nucleotide reactions: r1, r2, r1bac, r2bac. With bacterial catalysis
% Temperature is expressed in Celsius, time in days
% Reparametrized version: parameters are optimized in log-space (z = log(p))
clear variables
clc

% Load and organize experimental data
cd ../../../experimentalData
iComb   = 2;
purpose = 'Est';
expData = readExperimentalData(iComb,purpose);
cd ../knowledgeDrivenModels/vbmc/model07

% Parallel pool (process-based)
nExp   = length(expData);
pp = gcp('nocreate');
if isempty(pp)
    parpool('local', nExp); 
end

% Known parameters
T0Sh = -9.0;     % Temperature at which there is no growth [Celsius]. From Safety and quality issues in fish processing book Table 12.4. Pag 209
T0Ph = -9.0;     % Temperature at which there is no growth [Celsius]. From Safety and quality issues in fish processing book Table 12.4. Pag 209
A1   = 1.15e20;  % Pre-exponential coefficient in Arrhenius expression (reaction r1) [umol/(d·g-fish)]
A2   = 7.63e15;  % Pre-exponential coefficient in Arrhenius expression (reaction r2) [umol/(d·g-fish)]
R    = 8.314;    % Universal gas constant [J/(mol·K)]
kP   = [T0Sh, T0Ph, A1, A2, R];

% Bundle data for wrappers
otherData.expData = expData;
otherData.kP      = kP;
otherData.nu      = 4;         % Degrees of freedom for Student's t likelihood

% Load optimal value obtained in BADS
load('PEresults_model07.mat','Pbest')

% Number of physical parameters
nStates       = 6;                         % Number of states
otherData.npp = length(Pbest) - nStates;

% Build prior bounds around p_est
eps_pos = 1e-14;
lb      = max(0.7*Pbest, eps_pos);
plb     = max(0.8*Pbest, eps_pos);
pub     = max(1.2*Pbest, eps_pos);
ub      = max(1.4*Pbest, eps_pos);
% Parameters in z-space
z_lb  = log(lb);
z_plb = log(plb);
z_pub = log(pub);
z_ub  = log(ub);
z_est = log(Pbest);

% VBMC call
opts             = vbmc('defaults');
opts.MaxFunEvals = 1.0e5;

% Encourage stability in later stages
opts.NoiseShaping = false; 

% Save statistics
opts.SaveStats = true;

% VBMC function
fun_vbmc = @(z) costWrapper_07_logJoint(z, otherData, z_est);

% Call
VP = vbmc(fun_vbmc, z_est, z_lb, z_ub, z_plb, z_pub, opts);

% Posterior sampling
Xs_z = vbmc_rnd(VP, 1.0e4);
Xs   = exp(Xs_z);  

% Computation of first moments and other statistical quantities
post_mean   = mean(Xs);
post_std    = std(Xs);
post_median = median(Xs);
post_cov    = cov(Xs);
post_iqr    = quantile(Xs, [0.25 0.75]);

% Save the results
save('PEresults_model07_CI.mat', 'VP', 'post_mean', 'post_std', 'post_median', 's');

