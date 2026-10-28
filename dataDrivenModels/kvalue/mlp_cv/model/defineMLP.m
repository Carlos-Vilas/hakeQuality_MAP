function layers = defineMLP(numFeatures, hiddenSizes, dropRates)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Function to create a MLP with variable number of layers
% Input arguments:
% - hiddenSizes:  vector with sizes of each hidden layer
% - dropRates:    vector with dropout rate per layer
% - numFeatures:  number of predictors
% Output argument:
% - layers:  network layer architecture for the MLP regression model

layers = [
    featureInputLayer(numFeatures, ...
                      "Normalization","none", ...
                      "Name","input")
];

% Build hidden layers dynamically
numLayers = numel(hiddenSizes);

for i = 1:numLayers
    layers = [
        layers
        fullyConnectedLayer(hiddenSizes(i), "Name", sprintf("fc%d",i))
        reluLayer("Name", sprintf("relu%d",i))
        dropoutLayer(dropRates(i), "Name", sprintf("drop%d",i))
    ];
end

% Output layer
layers = [
    layers
    fullyConnectedLayer(1, "Name","fc_out")
    regressionLayer("Name","regression")
];

end