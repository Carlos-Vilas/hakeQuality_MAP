function [MAE, R2] = evaluateModelPerform(net, X, Y)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Function to evaluate the MAE and R2 values for a given prediction given by
% a neural network model
% Input arguments:
% - net:  NN model
% - X:    Predictors values
% - Y:    Response variable values
% Output argument:
% - MAE:  Mean absolute error
% - R2:   Coefficient of determination

% Prediction provided by the model
YPred = predict(net, X);

% Computation of the MAE
MAE = mean(abs(YPred - Y));

% R2 (coefficient of determination)
SS_res = sum((Y - YPred).^2);
SS_tot = sum((Y - mean(Y)).^2);
R2     = 1.0 - SS_res / SS_tot;

% Range of the true target values
minQ = min(Y);
maxQ = max(Y);

% Print metrics
fprintf('\n============================================================\n');
fprintf('RESULTS\n');
fprintf('------------------------------------------------------------\n');
fprintf('MAE          = %.4f\n', MAE);
fprintf('R2           = %.4f\n', R2);
fprintf('Range QIM    = [%.3f  ,  %.3f]\n', minQ, maxQ);
fprintf('============================================================\n');


end