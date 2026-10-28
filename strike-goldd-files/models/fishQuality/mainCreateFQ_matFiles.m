% Script to calls function "create_matFile_fishQuality" that creates a mat
% file used by Strike-Goldd to study structural identifiability of fish
% quality models considered in the manuscript:
%   Mathematical modeling of the degradation of biochemical and subjective
%   quality indicators in Hake fillets stored under modified atmosphere
%   packaging. Submitted to Journal of Food Enginering  
clear variables
clc

% Choose the input arguments
modName  = 'M01';    % Model name. Valid options 'M01', 'M02',..., 'M10'   
paramSet = 'Partial'; % Parameter set. 'Complete' for the whole set of 
                       % unknown parameters, 'Partial' for a subset

% Call the function. No output arguments. The function creates a mat file
% that is stores in this folder
create_matFile_fishQuality(modName, paramSet)