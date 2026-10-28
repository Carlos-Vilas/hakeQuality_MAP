% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Main script to select the best NN network among those trained with
% different values of the hyperparameters
clear variables
clc

% Set random number generator
rng(1);

% Load the models
fileNames = dir('results/*.mat');

% Initialize variables
nFiles      = length(fileNames);
mean_MAE_Va = zeros(1,nFiles);
se_MAE_Va   = zeros(1,nFiles);  
nPar        = zeros(1,nFiles);

% Load the results with different NN
for ii = 1:nFiles
    % Path to the ii mat file
    filePath = sprintf('%s/%s',fileNames(ii).folder,fileNames(ii).name);

    % Load file
    load(filePath, 'MAE_Va', 'nParams')

    % Data for 1-SE rule
    mean_MAE_Va(ii) = mean(MAE_Va);
    std_MAE_Va      = std(MAE_Va);
    se_MAE_Va(ii)   = std_MAE_Va/sqrt(length(MAE_Va));
    nPar(ii)        = nParams;
end

% Assign a large mean MAE to the files with zero meanMAE
mean_MAE_Va(mean_MAE_Va==0) = 1.0e5;

% 1-SE rule
[minMAE, idx]  = min(mean_MAE_Va);                    % Minimum value of mean MAE
se_minMAE      = se_MAE_Va(idx);                      % Standard error of MAEs for the minimum mean MAE
idxCand        = find(mean_MAE_Va<=minMAE+se_minMAE); % Models that fulfill the 1-SE rule
nParCand       = nPar(idxCand);                       % Number of trainable parameters in the candidates
[~, idxMinPar] = min(nParCand);
bestNN         = fileNames(idxCand(idxMinPar)).name;


% Load final Net and Check again
fprintf('\n\n\n\n\n\n\n\n\n\n')
fprintf('Best net\n')
fprintf('%s\n\n',bestNN)

% Load the best network
fileName = sprintf('%s/%s',fileNames(1).folder,bestNN);
load(fileName,'net')

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
cd ../dataDrivenModels/nucleotide/mlp_cv

% Build synthetic data used to train and select the best model and separate
% the data in modeling (original and synthetic) and test
nRep    = 50; % Number of synthetic repetitions to generate per time point
[Xmod_or, Ymod_or, Xmod, Ymod, Xtest, Ytest] = buildData(expDataMod, expDataTest, nRep);

% Number of features 
numFeatures = size(Xmod,2);

% Values used for data normaliztion
mu    = mean(Xmod,1);       % mean of training set
sigma = std(Xmod,0,1);      % std of training set
sigma(sigma == 0) = 1;      % avoid division by zero

%  Data normalization
Xmod    = (Xmod - mu)./sigma;
Xmod_or = (Xmod_or - mu)./sigma;
Xtest   = (Xtest - mu)./sigma;

% Evaluate model performance on training + test
[MAE_mod, R2_mod]    = evaluateModelPerform(net, Xmod_or, Ymod_or);
[MAE_test, R2_test]  = evaluateModelPerform(net, Xtest, Ytest);

% Plot results
Ymod_pred    = predict(net, Xmod);
Ymod_or_pred = predict(net, Xmod_or);
Ytest_pred   = predict(net, Xtest);

Ymin   = min([Ymod(:);Ymod_pred(:);Ymod_or(:);Ymod_or_pred(:);Ytest(:);Ytest_pred(:)]);
Ymax   = max([Ymod(:);Ymod_pred(:);Ymod_or(:);Ymod_or_pred(:);Ytest(:);Ytest_pred(:)]);
Ymin_b = min([Ymod_or(:);Ymod_or_pred(:);Ytest(:);Ytest_pred(:)]);
Ymax_b = max([Ymod_or(:);Ymod_or_pred(:);Ytest(:);Ytest_pred(:)]);

figure
hold on
scatter(Ymod, Ymod_pred, 60, [0.7 0.7 0.7], 'filled');
scatter(Ymod_or, Ymod_or_pred, 60, 'b', 'filled');
scatter(Ytest, Ytest_pred, 60, 'r', 'filled');
plot([Ymin, Ymax], [Ymin, Ymax],'k')
xlabel('Real QIM')
ylabel('Predicted QIM')
hold off


% Remove project subfolders
rmpath(fullfile(pwd, 'preprocessing'));
rmpath(fullfile(pwd, 'model'));
