function K = rbfKernel(X, Xref, gamma)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Function to compute the nonlinear RBF kernel
% Input arguments:
% - X:     Predictors values
% - Y:     Response variable values
% - gamma: RBF kernel width parameter 
% Output argument:
% - K:  RBF kernek

% Euclidean distance among samples of X and Xref
D = pdist2(X, Xref, 'euclidean').^2;

% RBF kernel computation
K = exp(-gamma * D);

end
