% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Main script to train and validate the RF models
clear variables
clc

rng(1);

% Select the number of trees
nTrees = 100;

% Select the minimum leaf size
mLeaf  = 8;

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

%  Random Forest Regression (TreeBagger)
%  Repeated CV 
nComb = 50; % Number of combinations train/validation
% Fraction of data that will be used for training
trainFrac = 0.8;

% Total number of sampling times
N = size(Xmod,1);

% Initialize arrays to store solutions
MAE_Tr  = zeros(nComb, 1);
R2_Tr   = zeros(nComb, 1);
MAE_Va  = zeros(nComb, 1);
R2_Va   = zeros(nComb, 1);

% Repeated Cross-Validation combinations
for ii = 1:nComb
    % Random permutations of the total number of data to generate the train
    % and the validation sets
    idx = randperm(N);
    nTrain = round(trainFrac*N);
    % Train set
    idxTrain = idx(1:nTrain);
    XTrain   = Xmod(idxTrain,:);
    YTrain   = Ymod(idxTrain,:);
    % Validation set
    idxVal = idx(nTrain+1:end);
    XVal   = Xmod(idxVal,:);
    YVal   = Ymod(idxVal,:);

    % Normalization
    XTrain = (XTrain - mu)./sigma;  % normalize training data
    XVal   = (XVal   - mu)./sigma;  % normalize validation data

    % Train the random forest model
    RFmodel = TreeBagger(nTrees, XTrain, YTrain, 'Method','regression',...
                         'MinLeafSize', mLeaf, 'OOBPrediction','off', ...
                         'Surrogate','off');

    % Evaluate performance metrics on the train and validation sets
    [MAE_Tr(ii), R2_Tr(ii)] = evaluateMetrics(RFmodel, XTrain, YTrain);
    [MAE_Va(ii), R2_Va(ii)] = evaluateMetrics(RFmodel, XVal, YVal);
end


% Save the model and the mean MAE for the validation set
saveName =  sprintf('results/RFmae_%02d',mLeaf);
save(saveName,'MAE_Va','MAE_Tr','mLeaf')

% Remove project subfolders
rmpath(fullfile(pwd, 'preprocessing'));

