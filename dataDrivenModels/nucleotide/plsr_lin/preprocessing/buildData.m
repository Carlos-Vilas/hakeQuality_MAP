function [Xmod_or, Ymod_or, Xmod, Ymod, Xtest, Ytest] = buildData(expDataMod, expDataTest, nRep)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%   
% Function to call the program that generates the synthetica data and
% to separate the data in three groups: modelling (original data),
% modeling (synthetic data), and test data
% Input arguments 
% - expDataMod:  Experimental data in the modeling set
% - expDataMod:  Experimental data in the test set
% - nRep:        Number of synthetic repetitions to generate per time point
% Output arguments
% - Xmod_or: Original data in the modeling set (predictors)
% - Ymod_or: Original data in the modeling set (response variable)
% - Xmod:    Synthetic data in the modeling set (predictors)
% - Ymod:    Synthetic data in the modeling set (response variable)
% - Xtest:   Original data in the test set (predictors)
% - Ytest:   Original data in the test set (response variable)
% Combine the modeling and test experimental data    

% Combine the modeling and test experimental data
expData = [expDataMod, expDataTest];

% Indexes of the experiments in the modeling and test sets inside expData
vExpMod  = 1:length(expDataMod);
vExpTest = length(expDataMod)+1:length(expDataMod)+length(expDataTest);

% Separate states
for iExp = 1:length(expData)
    Sh{iExp}  = expData{iExp}.Sh;
    Ph{iExp}  = expData{iExp}.Ph;
    TVB{iExp} = expData{iExp}.TVB*100; % TVB-N in g/100 g-fish
    IMP{iExp} = expData{iExp}.IMP;
    Ino{iExp} = expData{iExp}.Ino;
    Hx{iExp}  = expData{iExp}.Hx;
    QIM{iExp} = expData{iExp}.QIM;
end

% Identify the maximum number of replicates in all experiments and all variables. 
maxC = 0;
for iExp = union(vExpMod,vExpTest)
    maxC_Sh  = size(Sh{iExp},2); 
    maxC_Ph  = size(Ph{iExp},2);
    maxC_TVB = size(TVB{iExp},2);
    maxC_IMP = size(IMP{iExp},2);
    maxC_Ino = size(Ino{iExp},2);
    maxC_Hx = size(Hx{iExp},2);
    maxC_QIM = size(QIM{iExp},2);
    maxC_all = max([maxC_Sh, maxC_Ph, maxC_TVB, maxC_IMP, maxC_Ino, maxC_Hx, maxC_QIM]);
    if maxC_all > maxC
        maxC = maxC_all;
    end
end
% Fill the missing values for experiments with lower number of replicates
% with NaN
for iExp = union(vExpMod,vExpTest)
    ncSh = size(Sh{iExp},2);
    if ncSh<maxC
        Sh{iExp}(:,ncSh+1:maxC) = NaN;
    end
    ncPh = size(Ph{iExp},2);
    if ncPh<maxC
        Ph{iExp}(:,ncPh+1:maxC) = NaN;
    end
    ncTVB = size(TVB{iExp},2);
    if ncTVB<maxC
        TVB{iExp}(:,ncTVB+1:maxC) = NaN;
    end
    ncIMP = size(IMP{iExp},2);
    if ncIMP<maxC
        IMP{iExp}(:,ncIMP+1:maxC) = NaN;
    end
    ncIno = size(Ino{iExp},2);
    if ncIno<maxC
        Ino{iExp}(:,ncIno+1:maxC) = NaN;
    end
    ncHx = size(Hx{iExp},2);
    if ncHx<maxC
        Hx{iExp}(:,ncHx+1:maxC) = NaN;
    end
    ncQIM = size(QIM{iExp},2);
    if ncQIM<maxC
        QIM{iExp}(:,ncQIM+1:maxC) = NaN;
    end
end

% Initialize output arrays
Shmod  = [];
Phmod  = [];
TVBmod = [];
IMPmod = [];
Inomod = [];
Hxmod = [];
QIMmod = [];

% Remove rows with all columns NaN in the modeling set
for iExp = vExpMod
    % Keep only the data at times where all variables are defined
    tSh  = expData{iExp}.tsSh;
    tPh  = expData{iExp}.tsPh;
    tTVB = expData{iExp}.tsTVB;
    tIMP = expData{iExp}.tsIMP;
    tIno = expData{iExp}.tsIno;
    tHx  = expData{iExp}.tsHx;
    tQIM = expData{iExp}.tsQIM;
    tAll = intersect(intersect(intersect(intersect(intersect(intersect(tSh, tPh), tTVB), tIMP), tIno),tHx),tQIM);
    % Find the locations of common time values in each variable
    [~, iSh, ~]  = intersect(tSh, tAll);
    [~, iPh, ~]  = intersect(tPh, tAll);
    [~, iTVB, ~] = intersect(tTVB, tAll);
    [~, iIMP, ~] = intersect(tIMP, tAll);
    [~, iIno, ~] = intersect(tIno, tAll);
    [~, iHx, ~]  = intersect(tHx, tAll);
    [~, iQIM, ~] = intersect(tQIM, tAll);
     
    % Keep only the values at the common locations
    Shmod  = [Shmod; Sh{iExp}(iSh,:)];
    Phmod  = [Phmod; Ph{iExp}(iPh,:)];
    TVBmod = [TVBmod; TVB{iExp}(iTVB,:)]; 
    IMPmod = [IMPmod; IMP{iExp}(iIMP,:)];
    Inomod = [Inomod; Ino{iExp}(iIno,:)];
    Hxmod = [Hxmod; Hx{iExp}(iHx,:)];
    QIMmod = [QIMmod; QIM{iExp}(iQIM,:)];

end

% Original data
Xmod_or = [mean(Shmod,2,'omitnan'), mean(Phmod,2,'omitnan'), mean(TVBmod,2,'omitnan'), mean(IMPmod,2,'omitnan'), mean(Inomod,2,'omitnan'), mean(Hxmod,2,'omitnan')];
Ymod_or = mean(QIMmod,2,'omitnan');

% Generate synthetic data for the modeling set
Shsynth  = generateSyntheticData(Shmod, nRep);
Phsynth  = generateSyntheticData(Phmod, nRep);
TVBsynth = generateSyntheticData(TVBmod, nRep);
IMPsynth = generateSyntheticData(IMPmod, nRep);
Inosynth = generateSyntheticData(Inomod, nRep);
Hxsynth  = generateSyntheticData(Hxmod, nRep);
QIMsynth = generateSyntheticData(QIMmod, nRep);

% Put together the explanatory variables
Xmod = [Shsynth(:), Phsynth(:), TVBsynth(:), IMPsynth(:), Inosynth(:), Hxsynth(:)];
Ymod = QIMsynth(:);

% Generate the test set (no synthetic data)
% Initialize output arrays
Shtest  = [];
Phtest  = [];
TVBtest = [];
IMPtest = [];
Inotest = [];
Hxtest = [];
QIMtest = [];

% Remove rows with all columns NaN in the test set
for iExp = vExpTest

    % Keep only the data at times where all variables are defined
    tSh  = expData{iExp}.tsSh;
    tPh  = expData{iExp}.tsPh;
    tTVB = expData{iExp}.tsTVB;
    tIMP = expData{iExp}.tsIMP;
    tIno = expData{iExp}.tsIno;
    tHx  = expData{iExp}.tsHx;
    tQIM = expData{iExp}.tsQIM;
    tAll = intersect(intersect(intersect(intersect(intersect(intersect(tSh, tPh), tTVB), tIMP), tIno),tHx),tQIM);
    % Find the locations of common time values in each variable
    [~, iSh, ~]  = intersect(tSh, tAll);
    [~, iPh, ~]  = intersect(tPh, tAll);
    [~, iTVB, ~] = intersect(tTVB, tAll);
    [~, iIMP, ~] = intersect(tIMP, tAll);
    [~, iIno, ~] = intersect(tIno, tAll);
    [~, iHx, ~]  = intersect(tHx, tAll);
    [~, iQIM, ~] = intersect(tQIM, tAll);

    % Keep only rows with at least one value different for NaN
    Shtest  = [Shtest; Sh{iExp}(iSh,:)];
    Phtest  = [Phtest; Ph{iExp}(iPh,:)];
    TVBtest = [TVBtest; TVB{iExp}(iTVB,:)]; 
    IMPtest = [IMPtest; IMP{iExp}(iIMP,:)];
    Inotest = [Inotest; Ino{iExp}(iIno,:)];
    Hxtest = [Hxtest; Hx{iExp}(iHx,:)];
    QIMtest = [QIMtest; QIM{iExp}(iQIM,:)];

end

Xtest = [mean(Shtest,2,'omitnan'), mean(Phtest,2,'omitnan'), mean(TVBtest,2,'omitnan'), mean(IMPtest,2,'omitnan'), mean(Inotest,2,'omitnan'), mean(Hxtest,2,'omitnan')];
Ytest = mean(QIMtest,2,'omitnan');

% End function
end