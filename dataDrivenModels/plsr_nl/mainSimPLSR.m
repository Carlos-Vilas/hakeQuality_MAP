% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%  
% Main script to simulate the selected nonlinear partial least squares model

clear variables 
clc

% Fix random seed for reproducibility
rng(1); 

% Add project subfolders to MATLAB path
addpath(fullfile(pwd, 'preprocessing'));
addpath(fullfile(pwd, 'model'));

% Load and organize experimental data
cd ../../../experimentalData
iComb   = 2;
% Experiments used for training/validation
purpose = 'Est';
expDataMod = readExperimentalData(iComb,purpose);
% Experiments used for test
purpose = 'Val';
expDataTest = readExperimentalData(iComb,purpose);
cd ../dataDrivenModels/kvalue/plsr_nl

% Build synthetic data used to train and select the best model and separate
% the data in modeling (original and synthetic) and test
nRep    = 50; % Number of synthetic repetitions to generate per time point
[Xmod_or, Ymod_or, Xmod, Ymod, Xtest, Ytest] = buildData(expDataMod, expDataTest, nRep);

% Normalization parameters
mu    = mean(Xmod,1);   % mean of training set
sigma = std(Xmod,0,1);  % std of training set
sigma(sigma == 0) = 1;  % avoid division by zero

% Normalization in Xmod, Xmod_or and Xtest
Xmod     = (Xmod  - mu) ./ sigma;    % Modeling set with synthetic data
Xmod_or  = (Xmod_or  - mu) ./ sigma; % Original nodeling set
Xtest    = (Xtest - mu) ./ sigma;    % Test set

% User parameters
nComb    = 50;                       % Number of random repetitions for cross validation
maxComp  = min(40, size(Xmod,1)-1); % Maximum number of PLS components to test
N        = size(Xmod,1);

% Load the model
load('PLRSmodel_gamma_13.mat','beta','optComp','optGamma')

% Compute the Kernels for all sets
Kmod     = rbfKernel(Xmod, Xmod, optGamma);
Kmod_or  = rbfKernel(Xmod_or, Xmod, optGamma);
Ktest    = rbfKernel(Xtest, Xmod, optGamma);

% Evaluate the preformance of the model
[MAEmod, R2mod, Ymod_pred]          = evaluateMetrics(Xmod, Ymod, Kmod, beta);
[MAEmod_or, R2mod_or, Ymod_or_pred] = evaluateMetrics(Xmod_or, Ymod_or, Kmod_or, beta);
[MAEtest, R2test, Ytest_pred]       = evaluateMetrics(Xtest, Ytest, Ktest, beta);

% Print metrics
fprintf('\n============================================================\n');
fprintf('Performance metrics with the original model set\n');
fprintf('MAE          = %.4f\n', MAEmod_or);
fprintf('R2           = %.4f\n', R2mod_or);
fprintf('\n============================================================\n');
fprintf('Performance metrics with the test set\n');
fprintf('MAE          = %.4f\n', MAEtest);
fprintf('R2           = %.4f\n', R2test);

% Plot the results
% Ranges
YY = [Ymod_pred; Ymod; Ymod_or; Ymod_or_pred; Ytest;Ytest_pred];
Ymin = min(YY(:));
Ymax = max(YY(:));
% Create the figure
figure
hold on
scatter(Ymod, Ymod_pred, 60, [0.7 0.7 0.7], 'filled');
scatter(Ymod_or, Ymod_or_pred, 60, 'b', 'filled');
scatter(Ytest, Ytest_pred, 60, 'r', 'filled');
plot([Ymin, Ymax], [Ymin Ymax],'k')
xlabel('QIM measured')
ylabel('QIM predicted')
hold off

% Remove project subfolders
rmpath(fullfile(pwd, 'preprocessing'));
rmpath(fullfile(pwd, 'model'));
