% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging 
% Submitted to Journal of Food Engineering
%
% Main file to simulate model 01
clear variables
close all
clc

% Load and organize experimental data
cd ../../../experimentalData
iComb   = 2;
purpose = 'Est';
expData = readExperimentalData(iComb,purpose);
cd ../knowledgeDrivenModels/bads/model01

% Number of experiments
nExp   = length(expData);

% Known parameters
T0Sh = -9.0;     % Temperature at which there is no growth [Celsius]. From Safety and quality issues in fish processing book Table 12.4. Pag 209
T0Ph = -9.0;     % Temperature at which there is no growth [Celsius]. From Safety and quality issues in fish processing book Table 12.4. Pag 209
ATVB = 1.82e17;  % From: Yao, et al. 2011. Establishment of kinetic models based on electrical conductivity and freshness indicators for the forecasting of crucian carp (Carassius carassius) freshness
A1   = 1.15e20;  % Pre-exponential coefficient in Arrhenius expression (reaction r1) [umol/(d·g-fish)]
A2   = 7.63e15;  % Pre-exponential coefficient in Arrhenius expression (reaction r2) [umol/(d·g-fish)]
A3   = 2.29e18;  % Pre-exponential coefficient in Arrhenius expression (reaction r3) [umol/(d·g-fish)]
A4   = 4.84e12;  % Pre-exponential coefficient in Arrhenius expression (reaction r4) [umol/(d·g-fish)]
R    = 8.314;    % Universal gas constant [J/(mol·K)]
kP   = [T0Sh, T0Ph, ATVB, A1, A2, A3, A4, R];

% Load optimal values
load('PEresults_model01','Pbest','Jbest')

% Preallocation
tSim = cell(1, nExp);
ySim = cell(1, nExp);

% Model output indices
iSh = 1;  iPh = 2;  iTVB = 3;  iIMP = 4;  iIno = 5;  iHx = 6;

for iExp = 1:nExp

    % Simulate experiment iExp
    [tSim{iExp}, ySim{iExp}] = simulationModel01(Pbest(1:end-6), kP,...
        expData{iExp}.tT, expData{iExp}.T, expData{iExp}.sT,...
        iExp, nExp);

    % Plot simulation vs experimental data
    % Sh
    figure
    subplot(2,3,1)
    hold on
    errorbar(expData{iExp}.tsSh,mean(expData{iExp}.Sh,2,'omitnan'),std(expData{iExp}.Sh,0,2,'omitnan'),'o')
    plot(tSim{iExp}, ySim{iExp}(:,iSh),'-')
    hold off
    % Ph
    subplot(2,3,2)
    hold on
    errorbar(expData{iExp}.tsPh,mean(expData{iExp}.Ph,2,'omitnan'),std(expData{iExp}.Ph,0,2,'omitnan'),'o')
    plot(tSim{iExp}, ySim{iExp}(:,iPh),'-')
    hold off
    % TVB-N
    subplot(2,3,3)
    hold on
    errorbar(expData{iExp}.tsTVB,mean(expData{iExp}.TVB,2,'omitnan'),std(expData{iExp}.TVB,0,2,'omitnan'),'o')
    plot(tSim{iExp}, ySim{iExp}(:,iTVB),'-')
    hold off
    % IMP
    subplot(2,3,4)
    hold on
    errorbar(expData{iExp}.tsIMP,mean(expData{iExp}.IMP,2,'omitnan'),std(expData{iExp}.IMP,0,2,'omitnan'),'o')
    plot(tSim{iExp}, ySim{iExp}(:,iIMP),'-')
    hold off
    % Ino
    subplot(2,3,5)
    hold on
    errorbar(expData{iExp}.tsIno,mean(expData{iExp}.Ino,2,'omitnan'),std(expData{iExp}.Ino,0,2,'omitnan'),'o')
    plot(tSim{iExp}, ySim{iExp}(:,iIno),'-')
    hold off
    % Hpx
    subplot(2,3,6)
    hold on
    errorbar(expData{iExp}.tsHx,mean(expData{iExp}.Hx,2,'omitnan'),std(expData{iExp}.Hx,0,2,'omitnan'),'o')
    plot(tSim{iExp}, ySim{iExp}(:,iHx),'-')
    hold off

end


