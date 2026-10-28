% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Main script to train an MLP using trainNetwork
% One experiment is left unseen (test set)

clear variables

% Add project subfolders to MATLAB path
addpath(fullfile(pwd, 'preprocessing'));
addpath(fullfile(pwd, 'model'));

% Read hyperparameters from environment variables
nL = 3;   % Number of layers
nN = 6;   % Number of neurons per layer
wD = 0.0; % Weight decay
dV = 0.0; % Drop value

% Fixed hyperparameters
valPatience = 5.0e2;   % Early stopping patience
nEpochs     = 2.5e4;   % Maximum number of epochs
batchSize   = 50;      % Mini-batch size

% Hidden layers and dropout configuration
hiddenSizes = nN*ones(1, nL);
dropRates   = dV*ones(1, nL);

% Load and organize experimental data
cd ../../../experimentalData
iComb   = 2;
% Experiments used for training/validation
purpose = 'Est';
expDataMod = readExperimentalData(iComb,purpose);
% Experiments used for test
purpose = 'Val';
expDataTest = readExperimentalData(iComb,purpose);
cd ../dataDrivenModels/kvalue/mlp_cv

% Build synthetic data used to train and select the best model and separate
% the data in modeling (original and synthetic) and test
nRep    = 50; % Number of synthetic repetitions to generate per time point
[~, ~, Xmod, Ymod, ~, ~] = buildData(expDataMod, expDataTest, nRep);

% Values used for data normaliztion
mu    = mean(Xmod,1);       % mean of training set
sigma = std(Xmod,0,1);      % std of training set
sigma(sigma == 0) = 1;      % avoid division by zero

%  Data normalization (modeling set with synthetic values)
Xmod    = (Xmod - mu)./sigma;

% Internal validation split
nComb = 20; % Number of combinations train/validation

% Initialize arrays
MAE_Tr  = zeros(nComb, 1);
MAE_Va  = zeros(nComb, 1);

% Fraction of data that will be used for training
trainFrac = 0.8;       

% Number of features 
numFeatures = size(Xmod,2);
N           = size(Xmod,1);

% Define MLP using user-defined hyperparameters
layers = defineMLP(numFeatures, hiddenSizes, dropRates);

% Repeated Cross-Validation combinations
for ii = 1:nComb
    fprintf('Layers = %02d;\t Neurons = %02d;\t WD = %02f;\t DR = %02f;\t Comb = %02d of %02d\n', nL, nN, wD, dV, ii, nComb)
    % Random permutations of the total number of data to generate the train
    % and the validation sets
    idx = randperm(N);
    % Train set
    nTrain   = round(trainFrac*N);
    idxTrain = idx(1:nTrain);
    XTrain   = Xmod(idxTrain,:);
    YTrain   = Ymod(idxTrain,:);
    % Validation set
    idxVal = idx(nTrain+1:end);
    XVal   = Xmod(idxVal,:);
    YVal   = Ymod(idxVal,:);

    % Training options
    options = trainingOptions('adam', 'MaxEpochs', nEpochs, 'MiniBatchSize', batchSize, ...
                              'Shuffle', 'every-epoch', 'L2Regularization', wD, ...
                              'ValidationData', {XVal, YVal}, 'ValidationPatience', valPatience, ...
                              'Verbose', false, 'Plots', 'none');

    % Train network
    net = trainNetwork(XTrain, YTrain, layers, options);

    % Evaluate on training and validation
    [MAE_Tr(ii), ~] = evaluateModelPerform(net, XTrain, YTrain);
    [MAE_Va(ii), ~] = evaluateModelPerform(net, XVal, YVal);
end

% Number of model parameters (complexity)
nParams = 0;
for i = 1:numel(net.Layers)
    layer = net.Layers(i);

    % Fully connected layers
    if isprop(layer, 'Weights')
        nParams = nParams + numel(layer.Weights);
    end

    if isprop(layer, 'Bias')
        nParams = nParams + numel(layer.Bias);
    end

    % Batch normalization layers
    if isa(layer, 'nnet.cnn.layer.BatchNormalizationLayer')
        nParams = nParams + numel(layer.Scale); 
        nParams = nParams + numel(layer.Offset);
    end
end

% Train the model using the whole set of synthetic data 
options = trainingOptions('adam', 'MaxEpochs', nEpochs, 'MiniBatchSize', batchSize, ...
                          'Shuffle', 'every-epoch', 'L2Regularization', wD, ...
                          'Verbose', false, 'Plots', 'none');
net = trainNetwork(Xmod, Ymod, layers, options);

% Save the model
saveName = sprintf('results/finalMLPmodel_%02d_%02d_%02d_%02d',nL,nN,int32(wD*100),int32(dV*100));
save(saveName,'nParams','MAE_Va','MAE_Tr','net')


% Remove project subfolders
rmpath(fullfile(pwd, 'preprocessing'));
rmpath(fullfile(pwd, 'model'));
