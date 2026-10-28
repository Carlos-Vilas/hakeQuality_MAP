% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Main script to perform the impurity criteria analysis
clear variables
clc
rng(1);

% Select the optimal number of trees (computed using mainSelectRF.m)
nTrees = 100;

% Use the selected minimum leaf size (computed using mainSelectRF.m)
mLeaf = 9;

% Predictor names
predNames = {'Sh', 'Ps', 'TVB-N', '$K_I$-value'};

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
cd ../dataDrivenModels/kvalue/randomForest

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

% Initialize arrays for predictor importance analysis
nPred     = length(predNames);
rawImpCV  = zeros(nComb, nPred);
normImpCV = zeros(nComb, nPred);

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
    XTrain = (XTrain-mu)./sigma; % normalize training data
    XVal   = (XVal-mu)./sigma;   % normalize validation data

    % Train the random forest model
    RFmodel = TreeBagger(nTrees, XTrain, YTrain, 'Method', 'regression',...
                         'MinLeafSize', mLeaf, 'OOBPrediction', 'off', ...
                         'Surrogate', 'off');

    % Calculate impurity-based predictor importance
    rawImpCV(ii, :) = RFmodel.DeltaCriterionDecisionSplit;

    % Normalize the importance values within each fitted model
    normImpCV(ii, :) = rawImpCV(ii, :)/sum(rawImpCV(ii, :));
end

% Calculate importance statistics across CV repetitions
meanRawImp  = mean(rawImpCV, 1);
stdRawImp   = std(rawImpCV, 0, 1);
meanNormImp = mean(normImpCV, 1)*100;

% Sort predictors according to their mean importance
[meanImpSort, fRank] = sort(meanRawImp, 'descend');

% Print importance stability results
fprintf('Impurity-based predictor importance stability:\n\n')
fprintf('Predictor     Raw importance  Std. importance Norm. importance\n')
for ii = 1:nPred
    % Print information in screen
    fprintf('%-12s  %-14.2e  %-16.2e  %-10.2f \n', predNames{fRank(ii)},...
        meanRawImp(fRank(ii)), stdRawImp(fRank(ii)), meanNormImp(fRank(ii)));
end

% Plot mean importance and variability across CV repetitions
figure
%hold on
errorbar(1:nPred, meanRawImp(fRank), stdRawImp(fRank), 'k.', 'LineWidth', 1.2);
%hold off
xlabel('Predictor')
ylabel('Mean relative impurity-based importance')
xticks(1:nPred)
xticklabels(predNames(fRank))
grid off
box off
ax = gca;
ax.TickLabelInterpreter = 'latex';


% Remove project subfolders
rmpath(fullfile(pwd, 'preprocessing'));