% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% This file generates the C code "amigoRHS_model07.c" and compiles it
% together with the CVODEs library to create the MATLAB mex file to
% simulate Model 07. Such model considers:
% - Baranyi and Roberts model for bacterial growth
% - Power-Law model (only Ph, no effect of Sh) for TVB-N formation
% - Most important biochemical routes (r1, r2) for nucleotide degradation
%   WITH bacterial catalysis
clear variables
clear mex
warning off
clc

% Name of parameters
par_names{1}  = 'pbSh';
par_names{2}  = 'pT0Sh';
par_names{3}  = 'pyShmax';
par_names{4}  = 'pbPh';
par_names{5}  = 'pT0Ph';
par_names{6}  = 'pyPhmax';
par_names{7}  = 'pKPh';
par_names{8}  = 'pNPh';
par_names{9}  = 'pEa1';
par_names{10} = 'pA1';
par_names{11} = 'pEa2';
par_names{12} = 'pA2';
par_names{13} = 'pRR'  ;
par_names{14} = 'pk1bac'; 
par_names{15} = 'pk2bac' ;

% Name of the state variables
state_var{1} = 'ySh';
state_var{2} = 'yPh';
state_var{3} = 'yTVB';
state_var{4} = 'yIMP';
state_var{5} = 'yIno';
state_var{6} = 'yHpx';

% Name of inputs
stimulus_names{1} = 'sT';

% Number of variables
nstates   = length(state_var);
npar      = length(par_names);
nstimulus = length(stimulus_names);

% WRITING THE FORTRAN FUNCTIONS
fid  = fopen('C_src4Amigo/amigoRHS_model07.c','wt');

% Includes
fprintf(fid,'#include <amigoRHS.h>\n');
fprintf(fid,'\n');
fprintf(fid,'#include <math.h>\n');
fprintf(fid,'\n');
fprintf(fid,'#include <amigoJAC.h>\n');
fprintf(fid,'\n');
fprintf(fid,'#include <amigoSensRHS.h>\n');
fprintf(fid,'\n');
fprintf(fid,'#include <amigo_terminate.h>\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

% Definition of the state variables
fprintf(fid,'\n/* Definition of the state variables */\n');
for ii = 1:nstates
    fprintf(fid,'#define %s Ith(y,%u)\n',state_var{ii},ii-1);
end
fprintf(fid,'#define iexp amigo_model->exp_num\n'); % experiments

% Definition of the state variables derivatives
fprintf(fid,'\n/* Definition of the state variables derivatives*/\n');
for ii = 1:nstates
    fprintf(fid,'#define d%s Ith(ydot,%u)\n',state_var{ii},ii-1);
end

% Definition of parameters
fprintf(fid,'\n/* Definition of the model parameters*/\n');
for ii = 1:npar
    fprintf(fid,'#define %s (*amigo_model).pars[%u]\n',par_names{ii},ii-1);
end

% Definition of stimulus
fprintf(fid,'\n/* Definition of the stimulus*/\n');
for ii = 1:nstimulus
    fprintf(fid,'#define %s ((*amigo_model).controls_v[%u][(*amigo_model).index_t_stim]+(t-(*amigo_model).tlast)*(*amigo_model).slope[%u][(*amigo_model).index_t_stim])\n',stimulus_names{ii}, ii-1, ii-1);
end
fprintf(fid,'\n');

% Definition of algebraic variables
fprintf(fid,'\n\t/* Definition of algebraic variables */\n');
fprintf(fid,'\tdouble vmuSqrtSh;\n');
fprintf(fid,'\tdouble vmuSh;\n');
fprintf(fid,'\tdouble vSdiffSh;\n');
fprintf(fid,'\tdouble vmuSqrtPh;\n');
fprintf(fid,'\tdouble vmuPh;\n');
fprintf(fid,'\tdouble vSdiffPh;\n');
fprintf(fid,'\tdouble vShCFU;\n');
fprintf(fid,'\tdouble vPhCFU;\n');
fprintf(fid,'\tdouble vBacCFU;\n');
fprintf(fid,'\tdouble vk1;\n');
fprintf(fid,'\tdouble vr1;\n');
fprintf(fid,'\tdouble vk2;\n');
fprintf(fid,'\tdouble vr2;\n');
fprintf(fid,'\tdouble vr1Bac;\n');
fprintf(fid,'\tdouble vr2Bac;\n');
fprintf(fid,'\n');

% Write the RHS of the system
fprintf(fid,'\n/* Right hand side of the system */\n');
fprintf(fid,'int amigoRHS(realtype t, N_Vector y, N_Vector ydot, void *data){ \n');
fprintf(fid,'\tAMIGO_model* amigo_model=(AMIGO_model*)data;\n');
fprintf(fid,'\tctrlcCheckPoint(__FILE__, __LINE__);\n');

% Algebraic equations
fprintf(fid,'\n\t/* Algebraic equations */\n');
% Computation of intermediate variables for bacteria
% Sh
fprintf(fid,'\n\t/* Growth rate Sh*/\n');
fprintf(fid,'\tvmuSqrtSh = pbSh*(sT-pT0Sh);\n');
fprintf(fid,'\tvmuSh = pow(vmuSqrtSh,2.0);\n');
fprintf(fid,'\n\t/* Difference between the current bacteria (Sh) load and the maximum */\n');
fprintf(fid,'\tvSdiffSh = ySh-pyShmax;\n');
fprintf(fid,'\n\t/* Shewanella concentration in CFU */\n');
fprintf(fid,'\tvShCFU = pow(10.0,ySh);\n');
% Ph
fprintf(fid,'\n\t/* Growth rate Ph*/\n');
fprintf(fid,'\tvmuSqrtPh = pbPh*(sT-pT0Ph);\n');
fprintf(fid,'\tvmuPh = pow(vmuSqrtPh,2.0);\n');
fprintf(fid,'\n\t/* Difference between the current bacteria (Ph) load and the maximum */\n');
fprintf(fid,'\tvSdiffPh = yPh-pyPhmax;\n');
fprintf(fid,'\n\t/* Photobacterium concentration in CFU */\n');
fprintf(fid,'\tvPhCFU = pow(10.0,yPh);\n');
fprintf(fid,'\n');

% Total bacteria
fprintf(fid,'\n\t/* Total bacteria in CFU */\n');
fprintf(fid,'\tvBacCFU = vShCFU + vPhCFU;\n');
fprintf(fid,'\n');

% Computation of intermediate variables (reaction rates) for nucleotides
fprintf(fid,'\n\t/* Nucleotide reaction rates (enzymatic)*/\n');
fprintf(fid,'\tvk1 = pA1*exp(-(pEa1/(pRR*(sT+273.15))));\n');
fprintf(fid,'\tvr1 = vk1*yIMP;\n');
fprintf(fid,'\tvk2 = pA2*exp(-(pEa2/(pRR*(sT+273.15))));\n');
fprintf(fid,'\tvr2 = vk2*yIno;\n');
fprintf(fid,'\n\t/* Nucleotide reaction rates (bacteria)*/\n');
fprintf(fid,'\tvr1Bac = pk1bac*vBacCFU*yIMP;\n');
fprintf(fid,'\tvr2Bac = pk2bac*vBacCFU*yIno;\n');

% ODE system
fprintf(fid,'\n\t/* ODE system*/\n');
fprintf(fid,'\tdySh  = vmuSh*(1.0-pow(10.0,vSdiffSh));\n');
fprintf(fid,'\tdyPh  = vmuPh*(1.0-pow(10.0,vSdiffPh));\n');
fprintf(fid,'\tdyTVB = pKPh*pow(vPhCFU,pNPh);\n');
fprintf(fid,'\tdyIMP = -(vr1 + vr1Bac);\n');
fprintf(fid,'\tdyIno = vr1 + vr1Bac - (vr2 + vr2Bac);\n');
fprintf(fid,'\tdyHpx =  vr2 + vr2Bac;\n');
  
fprintf(fid,'\n');
fprintf(fid,'\treturn(0);\n');
fprintf(fid,'\n');
fprintf(fid,'}\n\n');

% Jacobian of the system
fprintf(fid,'\n/* Jacobian of the system*/\n');
fprintf(fid,'int amigoJAC(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){\n');
fprintf(fid,'\tAMIGO_model* amigo_model=(AMIGO_model*)user_data;\n');
fprintf(fid,'\tctrlcCheckPoint(__FILE__, __LINE__); \n');
fprintf(fid,'\n');
fprintf(fid,'\treturn(0);\n');
fprintf(fid,'\n');
fprintf(fid,'}\n\n');

% RHS sensitivity
fprintf(fid,'\n/* RHS sensitivity*/\n');
fprintf(fid,'int amigoSensRHS(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){ \n');
fprintf(fid,'\tAMIGO_model* amigo_model=(AMIGO_model*)data;\n');
fprintf(fid,'\n');
fprintf(fid,'\treturn(0);\n');
fprintf(fid,'\n');
fprintf(fid,'}\n\n');

fclose(fid);

% Mex compilation and generation
mex -compatibleArrayDims  -IC_src4Amigo/include/include_amigo -IC_src4Amigo/include/include_cvodes -lmwblas -lmwlapack ...
     -output mex_model07...
     C_src4Amigo/amigoRHS_model07.c ...
     C_src4Amigo/src/src_amigo/AMIGO_model.c ...
     C_src4Amigo/src/src_amigo/amigo_terminate.c ...
     C_src4Amigo/src/src_amigo/simulate_amigo_model.c ...
     C_src4Amigo/src/src_interface/interface_with_matlab.c ...
     C_src4Amigo/src/src_cvodes/cvodea.c ...
     C_src4Amigo/src/src_cvodes/cvodea_io.c ...
     C_src4Amigo/src/src_cvodes/cvodes.c ...
     C_src4Amigo/src/src_cvodes/cvodes_band.c ...
     C_src4Amigo/src/src_cvodes/cvodes_bandpre.c ...
     C_src4Amigo/src/src_cvodes/cvodes_bbdpre.c ...
     C_src4Amigo/src/src_cvodes/cvodes_dense.c ...
     C_src4Amigo/src/src_cvodes/cvodes_diag.c ...
     C_src4Amigo/src/src_cvodes/cvodes_direct.c ...
     C_src4Amigo/src/src_cvodes/cvodes_io.c ...
     C_src4Amigo/src/src_cvodes/cvodes_lapack.c ...
     C_src4Amigo/src/src_cvodes/cvodes_spbcgs.c ...
     C_src4Amigo/src/src_cvodes/cvodes_spgmr.c ...
     C_src4Amigo/src/src_cvodes/cvodes_spils.c ...
     C_src4Amigo/src/src_cvodes/cvodes_sptfqmr.c ...
     C_src4Amigo/src/src_cvodes/fnvector_serial.c ...
     C_src4Amigo/src/src_cvodes/nvector_serial.c ...
     C_src4Amigo/src/src_cvodes/sundials_band.c ...
     C_src4Amigo/src/src_cvodes/sundials_dense.c ...
     C_src4Amigo/src/src_cvodes/sundials_direct.c ...
     C_src4Amigo/src/src_cvodes/sundials_iterative.c ...
     C_src4Amigo/src/src_cvodes/sundials_math.c ...
     C_src4Amigo/src/src_cvodes/sundials_nvector.c ...
     C_src4Amigo/src/src_cvodes/sundials_spbcgs.c ...
     C_src4Amigo/src/src_cvodes/sundials_spgmr.c ...
     C_src4Amigo/src/src_cvodes/sundials_sptfqmr.c ;
