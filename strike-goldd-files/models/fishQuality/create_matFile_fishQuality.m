function create_matFile_fishQuality(modName, paramSet)
% This file creates the mat file used by Strike-goldd to analyze the
% structural identifiability of the fish quality models considered in the
% manuscript:
%   Mathematical modeling of the degradation of biochemical and subjective
%   quality indicators in Hake fillets stored under modified atmosphere
%   packaging. Submitted to Journal of Food Enginering
%
% Input arguments:
%  - modName: string that indicates the name of the model, i.e., 'M01',
%             'M02',..., 'M10'
%  - paramSet: string that indicates whether to consider all the unknown
%              model parameters to study structural identifiability of a
%              subset. To consider all parameters, choose 'Complete',
%              otherwise, choose 'Partial'. In the former some of the
%              parameters are structurally unidentifiable
%
% There are no output arguments. The function creates a mat file

%% Part of the code common to all models

% State variables
syms ySh yPh yTVB yIMP yIno yHpx
x = [ySh; yPh; yTVB; yIMP; yIno; yHpx];

% Measured variables (all states are measured)
h = [ySh; yPh; yTVB; yIMP; yIno; yHpx];

% One input (temperature)
syms sT;

% Parameters common to all models
syms pbSh pT0Sh pyShmax pbPh pT0Ph pyPhmax pEa1 pA1 pEa2 pA2 pEa3 pA3 pEa4 pA4 pRR pDD

% Auxiliary functions (growth rates, reaction rates, etc) common to all
% models
% Bacterial growth rates
vmuSh  = pbSh^2.0*(sT-pT0Sh);
vmuPh  = pbPh^2.0*(sT-pT0Ph);
% Bacterial concentration in CFU/g-fish
vShCFU  = 10.0^ySh;
vPhCFU  = 10.0^yPh;
vBacCFU = vShCFU + vPhCFU;
% Arrhenius equation in nucleotide degradation
vk1    = pA1*exp(-(pEa1/(pRR*(sT+273.15))));
vk2    = pA2*exp(-(pEa2/(pRR*(sT+273.15))));
vk3    = pA3*exp(-(pEa3/(pRR*(sT+273.15))));
vk4    = pA4*exp(-(pEa4/(pRR*(sT+273.15))));
% Reaction rates (first-order) in nucleotide degradation
vr1    = vk1*yIMP;
vr2    = vk2*yIno;
vr3    = vk3*yIMP;
vr4    = vk4*yIMP;

% Baranyi and Roberts equations: used in all models
modelSh = vmuSh*(1.0-10.0^(ySh-pyShmax));
modelPh = vmuPh*(1.0-10.0^(yPh-pyPhmax));

% Initial conditions
ics       = [];
known_ics = [];


%% Unknown parameters in the different models
switch modName
    case{'M01'}
        % Model 01 which considers:
        % - Baranyi and Roberts model for bacterial growth
        % - First-order model (without bacterial influence) for TVB-N formation
        % - Complete biochemical route for nucleotide degradation WITHOUT bacterial catalysis
        % Parameters particular of this model
        syms pATVB pEaTVB
        switch paramSet
            case{'Complete'}
                % 17 unknown
                p = [pbSh; pT0Sh; pyShmax; pbPh; pT0Ph; pyPhmax; pATVB; pEaTVB;...
                    pEa1; pA1; pEa2; pA2; pEa3; pA3; pEa4; pA4; pDD];
            case{'Partial'}
                % 10 unknown
                p = [pbSh; pyShmax; pbPh; pyPhmax; pEaTVB; pEa1; pEa2; pEa3; pEa4; pDD];
            otherwise
                error('MyFunction:InvalidParamSet', ...
                    'Invalid paramSet "%s". Valid options are: ''Complete'', ''Partial'' ', ...
                    paramSet);
        end

    case{'M02'}
        % Model 02 which considers:
        % - Baranyi and Roberts model for bacterial growth
        % - Power-law model (both bacteria: Sh, Ph) for TVB-N formation
        % - Complete biochemical route for nucleotide degradation WITHOUT bacterial catalysis
        % Parameters particular of this model
        syms pKSh pNSh pKPh pNPh
        switch paramSet
            case{'Complete'}
                % 19 unknown
                p = [pbSh; pT0Sh; pyShmax; pbPh; pT0Ph; pyPhmax; pKSh; pNSh;...
                    pKPh; pNPh; pEa1; pA1; pEa2; pA2; pEa3; pA3; pEa4; pA4; pDD];
            case{'Partial'}
                % 13 unknown
                p = [pbSh; pyShmax; pbPh; pyPhmax; pKSh; pNSh; pKPh; pNPh;...
                    pEa1; pEa2; pEa3; pEa4; pDD];
            otherwise
                error('MyFunction:InvalidParamSet', ...
                    'Invalid paramSet "%s". Valid options are: ''Complete'', ''Partial'' ', ...
                    paramSet);
        end

    case{'M03'}
        % Model 03 which considers:
        % - Baranyi and Roberts model for bacterial growth
        % - Luedeking-Piret model (both bacteria: Sh, Ph) for TVB-N formation
        % - Complete biochemical route for nucleotide degradation WITHOUT bacterial catalysis
        % Parameters particular of this model
        syms pYSh pmSh pYPh pmPh
        switch paramSet
            case{'Complete'}
                % 19 unknown
                p = [pbSh; pT0Sh; pyShmax; pbPh; pT0Ph; pyPhmax; pYSh; pmSh;...
                    pYPh; pmPh; pEa1; pA1; pEa2; pA2; pEa3; pA3; pEa4; pA4; pDD];
            case{'Partial'}
                % 13 unknown
                p = [pbSh; pyShmax; pbPh; pyPhmax; pYSh; pmSh; pYPh; pmPh;...
                    pEa1; pEa2; pEa3; pEa4; pDD];
            otherwise
                error('MyFunction:InvalidParamSet', ...
                    'Invalid paramSet "%s". Valid options are: ''Complete'', ''Partial'' ', ...
                    paramSet);
        end
    case{'M04'}
        % Model 04 which considers:
        % - Baranyi and Roberts model for bacterial growth
        % - First-order model (without bacterial influence) for TVB-N formation
        % - Complete biochemical route for nucleotide degradation WITH bacterial catalysis
        % Parameters particular of this model
        syms pATVB pEaTVB pk1bac pk2bac
        switch paramSet
            case{'Complete'}
                % 19 unknown
                p = [pbSh; pT0Sh; pyShmax; pbPh; pT0Ph; pyPhmax; pATVB; pEaTVB;...
                    pEa1; pA1; pEa2; pA2; pEa3; pA3; pEa4; pA4; pk1bac; pk2bac; pDD];
            case{'Partial'}
                % 12 unknown
                p = [pbSh; pyShmax; pbPh; pyPhmax; pEaTVB; pEa1; pEa2; pEa3; pEa4; pk1bac; pk2bac; pDD];
            otherwise
                error('MyFunction:InvalidParamSet', ...
                    'Invalid paramSet "%s". Valid options are: ''Complete'', ''Partial'' ', ...
                    paramSet);
        end

    case{'M05'}
        % Model 05 which considers:
        % - Baranyi and Roberts model for bacterial growth
        % - Power-law model (both bacteria: Sh, Ph) for TVB-N formation
        % - Complete biochemical route for nucleotide degradation WITH bacterial catalysis
        % Parameters particular of this model
        syms pKSh pNSh pKPh pNPh pk1bac pk2bac
        switch paramSet
            case{'Complete'}
                % 21 unknown
                p = [pbSh; pT0Sh; pyShmax; pbPh; pT0Ph; pyPhmax; pKSh; pNSh;...
                    pKPh; pNPh; pEa1; pA1; pEa2; pA2; pEa3; pA3; pEa4; pA4; pk1bac; pk2bac; pDD];
            case{'Partial'}
                % 15 unknown
                p = [pbSh; pyShmax; pbPh; pyPhmax; pKSh; pNSh; pKPh; pNPh;...
                    pEa1; pEa2; pEa3; pEa4; pk1bac; pk2bac; pDD];
            otherwise
                error('MyFunction:InvalidParamSet', ...
                    'Invalid paramSet "%s". Valid options are: ''Complete'', ''Partial'' ', ...
                    paramSet);
        end

    case{'M06'}
        % Model 06 which considers:
        % - Baranyi and Roberts model for bacterial growth
        % - Luedeking-Piret model (both bacteria: Sh, Ph) for TVB-N formation
        % - Complete biochemical route for nucleotide degradation WITH bacterial catalysis
        % Parameters particular of this model
        syms pYSh pmSh pYPh pmPh pk1bac pk2bac
        switch paramSet
            case{'Complete'}
                % 21 unknown
                p = [pbSh; pT0Sh; pyShmax; pbPh; pT0Ph; pyPhmax; pYSh; pmSh;...
                    pYPh; pmPh; pEa1; pA1; pEa2; pA2; pEa3; pA3; pEa4; pA4;...
                    pk1bac; pk2bac; pDD];
            case{'Partial'}
                % 15 unknown
                p = [pbSh; pyShmax; pbPh; pyPhmax; pYSh; pmSh; pYPh; pmPh;...
                    pEa1; pEa2; pEa3; pEa4; pk1bac; pk2bac; pDD];
            otherwise
                error('MyFunction:InvalidParamSet', ...
                    'Invalid paramSet "%s". Valid options are: ''Complete'', ''Partial'' ', ...
                    paramSet);
        end

    otherwise
        error('MyFunction:InvalidModName', ...
            'Invalid modName "%s". Valid options are: ''M01'', ''M02'',...,''M10'' ', ...
            modName);
end

%% TVB-N models
switch modName
    case{'M01','M04'}
        % Arrhenius equation in TVB model
        kTVB   = pATVB*exp(-(pEaTVB/(pRR*(sT+273.15))));
        % Model TVB-N
        modelTVB = kTVB*yTVB;
    case{'M02','M05'}
        % Model TVB-N
        modelTVB = pKSh*vShCFU^pNSh + pKPh*vPhCFU^pNPh;
    case{'M03','M06'}
        % Model TVB-N
        modelTVB = pYSh*modelSh + pmSh*ySh + pYPh*modelPh + pmPh*yPh;
    otherwise
        error('MyFunction:InvalidModName', ...
            'Invalid modName "%s". Valid options are: ''M01'', ''M02'',...,''M10'' ', ...
            modName);
end



%% Nucleotide models
switch modName
    case{'M01','M02','M03'} % No bacteria catalysis in nucleotide model
        % Nucleotide
        modelIMP = -(vr1 + vr3 + vr4) - pDD*yIMP;
        modelIno = vr1 - vr2 - pDD*yIno;
        modelHpx = vr2 + vr4 - pDD*yHpx;
    case{'M04','M05','M06'} % Bacteria catalysis in nucleotide model
        % Bacterial catalysis reaction rates
        vr1bac = pk1bac*vBacCFU*yIMP;
        vr2bac = pk2bac*vBacCFU*yIno;
        % Nucleotide
        modelIMP = -(vr1 + vr1bac + vr3 + vr4) - pDD*yIMP;
        modelIno = vr1 + vr1bac - vr2 - vr2bac - pDD*yIno;
        modelHpx = vr2 + vr2bac + vr4 - pDD*yHpx;
    otherwise
        error('MyFunction:InvalidModName', ...
            'Invalid modName "%s". Valid options are: ''M01'', ''M02'',...,''M10'' ', ...
            modName);
end


% Dynamic equations
f = [modelSh;
    modelPh;
    modelTVB;
    modelIMP;
    modelIno;
    modelHpx];

% Save the matfile to be used in Strike-goldd
saveName = sprintf('%s_%s_structIdent',modName,paramSet);
save(saveName,'x','p','h','f','ics','known_ics');