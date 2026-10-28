% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%  
% Main script to train and validate the nonlinear partial least squares models
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

% Decay parameter in the RBF kernel
gamma0  = 1/size(Xmod,2);   % Default value
gamma_v = gamma0*logspace(-2,2,50);
gamma_v = union(gamma_v,gamma0);

% Normalization parameters
mu    = mean(Xmod,1);   % mean of training set
sigma = std(Xmod,0,1);  % std of training set
sigma(sigma == 0) = 1;  % avoid division by zero

% Normalization in Xmod, Xmod_or and Xtest
Xmod     = (Xmod  - mu)./sigma;    % Modeling set with synthetic data
Xmod_or  = (Xmod_or  - mu)./sigma; % Original nodeling set
Xtest    = (Xtest - mu)./sigma;    % Test set

% User parameters
nComb    = 50;                       % Number of random repetitions for cross validation
maxComp  = min(40, size(Xmod,1)-1); % Maximum number of PLS components to test
N        = size(Xmod,1);

% Fraction of data that will be used for training
trainFrac = 0.8;

% Initialize arrays
MAE_Tr = zeros(nComb,length(gamma_v),maxComp);
MAE_Va = zeros(nComb,length(gamma_v),maxComp);

% Create parallel pool
nProc = 4;
% Get existing parallel pool (if any)
pool = gcp('nocreate');
if isempty(pool)
    % No pool exists -> create a new one
    pool = parpool('local', nProc);
else
    % Pool exists -> check number of workers
    if pool.NumWorkers ~= nProc
        % Delete existing pool and create a new one
        delete(pool);
        pool = parpool('local', nProc);
    end
end

% Cross-validation combinations
for ii = 1:nComb
    fprintf('Combination %02d of %02d\n',ii,nComb)
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

    for jj = 1:length(gamma_v)

        % Build the RBF kernel for the train and validation sets
        KTrain = rbfKernel(XTrain, XTrain, gamma_v(jj));
        KVal   = rbfKernel(XVal,   XTrain, gamma_v(jj));   % always vs training set

        % Loop over number of PLS components
        parfor kk = 1:maxComp

            % Fit PLS model
            [XL, YL, XS, YS, beta, PCTVAR] = plsregress(KTrain, YTrain, kk);

            % Compute MAE
            [MAE_Va(ii,jj,kk), ~, ~]  = evaluateMetrics(XVal, YVal, KVal, beta);
        end
    end
end

% Select optimal value of the hyperparameters according to the 1-SE rule
% Compute mean and standard error across CV repetitions
mean_MAE_Va = squeeze(mean(MAE_Va, 1));   
std_MAE_Va  = squeeze(std(MAE_Va, 0, 1)); 
se_MAE_Va   = std_MAE_Va / sqrt(nComb);

% Find global minimum of mean MAE
[minMAE, idxMin] = min(mean_MAE_Va(:));
[igOpt, jcOpt]   = ind2sub(size(mean_MAE_Va), idxMin);

% Compute 1-SE threshold
MAE_threshold = minMAE + se_MAE_Va(igOpt, jcOpt);

% Apply hierarchical 1-SE rule
flag = 0;  % Flag to stop the search. 0 continue the search, 1 stop it
for jc = 1:size(mean_MAE_Va,2)        % loop over number of components
    for ig = 1:size(mean_MAE_Va,1)    % loop over gamma values
        if mean_MAE_Va(ig,jc) <= MAE_threshold
            idxGamma = ig;
            optComp  = jc;
            flag = 1;
            break
        end
    end
    if flag == 1
        break
    end
end
optGamma = gamma_v(idxGamma);

% Print information in screen
fprintf('\n\n\n\n\nOptimal number of components = %d\n', optComp);
fprintf('Optimal gamma = %d\n\n\n\n\n', optGamma);

% Compute the Kernels for all sets
Kmod     = rbfKernel(Xmod, Xmod, optGamma);
Kmod_or  = rbfKernel(Xmod_or, Xmod, optGamma);
Ktest    = rbfKernel(Xtest, Xmod, optGamma);

% Train the final model with the optimal number of parameters
[XL, YL, XS, YS, beta, PCTVAR] = plsregress(Kmod, Ymod, optComp);

% Save the model
saveName = sprintf('PLRSmodel_gamma_%02d',optComp);
save(saveName,'optComp','optGamma','beta','PCTVAR')

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
