function [MAE, R2] = evaluateMetrics(RFmodel, X, Y)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Function to compute the metrics for the evaluation of the model
% predictive capabilities
% Input arguments:
% - RFmodel: Random forest model
% - X:       Predictors values
% - Y:       Response variable values
% Output argument:
% - MAE:    Mean absolute error
% - R2:     Coefficient of determination

% Model predictions
YPred = predict(RFmodel, X);

% Metrics
MAE       = mean(abs(YPred - Y));
SS_res    = sum((Y - YPred).^2);
SS_tot    = sum((Y - mean(Y)).^2);
R2        = 1 - SS_res / SS_tot;

end

