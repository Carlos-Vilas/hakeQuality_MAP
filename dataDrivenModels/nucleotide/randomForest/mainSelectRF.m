% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Main script to select the best RF model
clear variables 
clc

rng(1);

% Select the best model
fileNames = dir('results/*.mat');
% Initialize variables
nFiles      = length(fileNames);
mean_MAE_Va = zeros(1,nFiles);
se_MAE_Va   = zeros(1,nFiles);
mLeafv      = zeros(1,nFiles);
for ii = 1:nFiles

    % Path to the ii mat file
    filePath = sprintf('%s/%s',fileNames(ii).folder,fileNames(ii).name);
    
    % Load the data
    load(filePath,'MAE_Va','mLeaf')

    % Data for 1-SE rule
    mean_MAE_Va(ii) = mean(MAE_Va);
    std_MAE_Va      = std(MAE_Va);
    se_MAE_Va(ii)   = std_MAE_Va/sqrt(length(MAE_Va));
    mLeafv(ii)      = mLeaf;
end

% 1-SE rule
[minMAE, idx]  = min(mean_MAE_Va);                    % Minimum value of mean MAE
se_minMAE      = se_MAE_Va(idx);                      % Standard error of MAEs for the minimum mean MAE
idxCand        = find(mean_MAE_Va<=minMAE+se_minMAE); % Models that fulfill the 1-SE rule
mLeafCand      = mLeafv(idxCand);                     % Number of trainable parameters in the candidates
[~, idxMinPar] = max(mLeafCand);
bestRF         = fileNames(idxCand(idxMinPar)).name;  % Best random forest

% Print information about the best model
fprintf('Best RF model:\n')
fprintf('Maximum number of leafs: %03d\n\n', mLeafv(idxCand(idxMinPar)))

% Select the number of trees
nTrees = 100;

% Minimum number of Leafs
mLeafs = mLeafv(idxCand(idxMinPar));

% Add project subfolders to MATLAB path
addpath(fullfile(pwd, 'preprocessing'));

% Load and organize experimental data
cd ../../../experimentalData
iComb   = 2;
% Experiments used for training/validation
purpose = 'Est';
expDataMod = readExperimentalData(iComb,purpose);
% Experiments used for test
purpose = 'Val';
expDataTest = readExperimentalData(iComb,purpose);
cd ../dataDrivenModels/nucleotide/randomForest

% Build synthetic data used to train and select the best model and separate
% the data in modeling (original and synthetic) and test
nRep    = 50; % Number of synthetic repetitions to generate per time point
[Xmod_or, Ymod_or, Xmod, Ymod, Xtest, Ytest] = buildData(expDataMod, expDataTest, nRep);

% Parameters for data normalization
mu    = mean(Xmod,1);       % mean of training set
sigma = std(Xmod,0,1);   % std of training set
sigma(sigma == 0) = 1;   % avoid division by zero

%  Normalization
Xmod    = (Xmod - mu)./sigma;
Xmod_or = (Xmod_or - mu)./sigma;
Xtest   = (Xtest - mu)./sigma;

% Train with the best set of hiperparameters and the synthetic data
RFmodel = TreeBagger(nTrees, Xmod, Ymod, 'Method','regression', ...
                     'MinLeafSize', mLeafs, 'OOBPrediction','off', ...
                     'Surrogate','off');
% Save the model
saveName = sprintf('bestRFmodel_%02dLeafs',mLeafs);
save(saveName,'RFmodel','mLeafs','mu','sigma')

% Evaluate the preformance of the model
[MAEmod, R2mod]       = evaluateMetrics(RFmodel, Xmod, Ymod);
[MAEmod_or, R2mod_or] = evaluateMetrics(RFmodel, Xmod_or, Ymod_or);
[MAEtest, R2test]     = evaluateMetrics(RFmodel, Xtest, Ytest);

% Print metrics
fprintf('\n============================================================\n');
fprintf('Performance metrics with the model set\n');
fprintf('MAE          = %.4f\n', MAEmod_or);
fprintf('R2           = %.4f\n', R2mod_or);
fprintf('\n============================================================\n');
fprintf('Performance metrics with the test set\n');
fprintf('MAE          = %.4f\n', MAEtest);
fprintf('R2           = %.4f\n', R2test);

% Plot the results
Ymod_pred    = predict(RFmodel, Xmod);
Ymod_or_pred = predict(RFmodel, Xmod_or);
Ytest_pred   = predict(RFmodel, Xtest);

% Ranges
YY = [Ymod_pred; Ymod; Ymod_or; Ymod_or_pred; Ytest;Ytest_pred];
Ymin = min(YY(:));
Ymax = max(YY(:));

% Plot
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