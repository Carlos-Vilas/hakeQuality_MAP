function [tT, ySim] = simulationModel03(eP, kP, tT, T, sT, iExp, nExp)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%    
% Function to simulate model 03
% Input arguments:
% - eP:   Parameters to be estimated
% - kP:   Known parameters. 
% - tT:   Measurement times for the control (temperature) [d]
% - T:    Values of the temperature at the measurement times [C]
% - sT:   Time-Temperature slopes [C/d]
% - ts:   Sampling times for the measurements of Photobacterium [d]
% - kk:   Relative number of experiment in this simulation
% Output arguments:
% - ySim: States at the integration times tT
%
% Note: Temperature measurement times (tT) will be used as simulation times
%       to avoid problems with too small intervals when using other times.

% Model known parameters 
T0Sh     = kP(1);
T0Ph     = kP(2);
A1       = kP(3);
A2       = kP(4);
A3       = kP(5);
A4       = kP(6);
R        = kP(7);

% Model unknown parameters
bSh   = eP(1); 
Shmax = eP(2);
bPh   = eP(3); 
Phmax = eP(4);
YSh   = eP(5); 
mSh   = eP(6); 
YPh   = eP(7); 
mPh   = eP(8); 
Ea1   = eP(9);
Ea2   = eP(10);
Ea3   = eP(11);
Ea4   = eP(12);
D     = eP(13);

% Total parameters
par     = [bSh, T0Sh, Shmax, bPh, T0Ph, Phmax, YSh, mSh, YPh, mPh,...
           Ea1, A1, Ea2, A2, Ea3, A3, Ea4, A4, R, D];
npar    = length(par);

% Initial conditions
nGP   = length(eP)-6*nExp; % Number of global parameters
ySh0  = eP(nGP+iExp);
yPh0  = eP(nGP+nExp+iExp);
yTVB0 = eP(nGP+2*nExp+iExp);
yIMP0 = eP(nGP+3*nExp+iExp);
yIno0 = eP(nGP+4*nExp+iExp);
yHpx0 = eP(nGP+5*nExp+iExp);
y0    = [ySh0, yPh0, yTVB0, yIMP0, yIno0, yHpx0];
nst   = length(y0);

% Number of simulation and control times, and number of controls
nt = length(tT);
nu = 1;

% IVP solver parameters:
rtol        = 1.0e-6;   % Relative integration tolerance
atol        = 1.0e-6;   % Absolute integration tolerance
jacobian    = 0;        % No use of the Jacobian 
Mstepsize   = 50;       % Maximum step size
Mnumsteps   = 1e7;      % Maximum number of steps
Metestfails = 5e3;      % Maximum error test fails
sens_an     = 0;        % No sensitivity analysis

% Call the mex file to numerically solve the model
[ySim, ~, ~] = mex_model03(npar, 0, par, ones(1,npar), tT(1), tT(nt), nt,...
                           tT, nst, y0, nu, nt, tT, T, sT, rtol, atol,...
                           Mstepsize, Mnumsteps, Metestfails, sens_an,...
                           jacobian, iExp-1, zeros(1,nst)); 


% End of the function
end
