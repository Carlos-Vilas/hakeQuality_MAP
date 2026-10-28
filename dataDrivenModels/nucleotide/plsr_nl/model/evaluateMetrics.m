function [MAE, R2, YPred] = evaluateMetrics(X, Y, K, beta)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Function to compute the metrics for the evaluation of the model
% predictive capabilities
% Input arguments:
% - X:    Predictors values
% - Y:    Response variable values
% - K:    RBF kernel
% - beta: Coefficients of the PLRS model
% Output argument:
% - MAE:    Mean absolute error
% - R2:     Coefficient of determination
% - Ypred:  Response variables values predicted by the model


% Prediction using the nonlinear PLRS model
YPred = [ones(size(X,1),1) K]*beta;

% Metrics
MAE       = mean(abs(YPred - Y));
SS_res    = sum((Y - YPred).^2);
SS_tot    = sum((Y - mean(Y)).^2);
R2        = 1 - SS_res / SS_tot;

end

