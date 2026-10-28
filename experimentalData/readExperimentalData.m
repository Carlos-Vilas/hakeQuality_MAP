function expData = readExperimentalData(iComb, purpose)
% Manuscript:
% Mathematical modeling of the degradation of biochemical and subjective
% quality indicators in Hake fillets stored under modified atmosphere packaging
% Submitted to Journal of Food Engineering
%
% Function to read, preprocess and arrange the experimental data used for
% model calibration and validation
%
% Input arguments:
% - iComb:   Cross-validation combination identifier. Four combinations are
%            available, each one using three experiments for parameter
%            estimation and the remaining experiment for validation
% - purpose: Data subset to retrieve:
%            'Est' for parameter estimation. 'Val' for model validation
%
% Output argument:
% - expData: Cell array containing the experimental data structures for the
%            selected experiments. Each structure includes:
%            Biochemical and microbiological measurements (Shewanella, Photobacterium, 
%            TVB-N, IMP, Ino, Hx and QIM)
%            Sampling times associated with each quality indicator
%            Temperature profile and corresponding time vector
%            Temperature slopes used as control inputs in the model

% Load the experimental data in the datase
conn = sqlite("hakeQuality_expData.db");

% Retrieve biological and chemical measurements
bioChemData = fetch(conn, "SELECT * FROM measurements");

% Retrieve temperature measurements
% datetime is stored as ISO 8601 TEXT
Tdata = fetch(conn, ...
    "SELECT * FROM temperature ORDER BY experiment_id, elapsed_days");
Tdata.datetime = datetime(Tdata.datetime, ...
    "InputFormat", "yyyy-MM-dd'T'HH:mm:ss");

close(conn)

% Separate variables and transform data to cell arrays
% Experiments available 
experiment_ids = unique(bioChemData.experiment_id);

% Biochemical variables
variable_ids = unique(bioChemData.variable_id);

% Determine largest experiment index
maxExperiment = 0;
for iExp = 1:numel(experiment_ids)
    experimentName   = char(experiment_ids(iExp));
    experimentNumber = sscanf(experimentName,'Experiment_%d');
    maxExperiment    = max(maxExperiment,experimentNumber);
end

% Pre-create one cell array per variable
for iVar = 1:numel(variable_ids)
    currentVariable = variable_ids(iVar);
    eval(sprintf('%s = cell(%d,1);',char(currentVariable),maxExperiment));
end

% Process experiment by experiment
for iExp = 1:numel(experiment_ids)
    % Extract the experiment number and the corresponding index
    experimentName   = experiment_ids(iExp);
    experimentNumber = sscanf(char(experimentName),'Experiment_%d');
    idxExperiment    = bioChemData.experiment_id == experimentName;

    % Data corresonding to experiment idxExperiment
    experimentDataAll = bioChemData(idxExperiment,:);

    % Sampling times for this experiment
    ts{experimentNumber} = unique(experimentDataAll.time_days);

    % Data for the measures variables
    for iVar = 1:numel(variable_ids)

        % Current variable name and data
        currentVariable = variable_ids(iVar);
        idxVariable     = experimentDataAll.variable_id == currentVariable;
        variableData    = experimentDataAll(idxVariable,:);

        % Unique times for this experiment and variable
        experimentTimes = unique(variableData.time_days);

        % Determine matrix width
        maxColumns = 0;
        for iTime = 1:numel(experimentTimes)
            idxTime    = variableData.time_days == experimentTimes(iTime);
            maxColumns = max(maxColumns,sum(idxTime));
        end

        % Create matrix
        M = nan(numel(experimentTimes),maxColumns);

        % Fill matrix
        for iTime = 1:numel(experimentTimes)
            idxTime = variableData.time_days == experimentTimes(iTime);
            valuesAtThisTime = variableData.value(idxTime);
            M(iTime,1:numel(valuesAtThisTime)) = valuesAtThisTime(:)';
        end

        % Store matrix in corresponding cell
        eval(sprintf('%s{%d} = M;', char(currentVariable), experimentNumber));

    end

    % Data for the temperature
    idxTemp    = strcmp(Tdata.experiment_id, experimentName);
    tcon_inter = Tdata.elapsed_days(idxTemp)';
    T_inter    = Tdata.temp_c(idxTemp)';

    % Remove extra values (those whose time is larger than ts{iExp}(end))
    idxCut = find(tcon_inter > ts{experimentNumber}(end));
    tcon_inter(idxCut) = [];
    T_inter(idxCut)    = [];

    % Ensure that the last sampling time is included in the temperature data
    if ts{experimentNumber}(end) ~= tcon_inter(end)
        tcon_inter(end) = ts{experimentNumber}(end);
    end

    % New sampling times for the control variables
    tnew = tcon_inter;

    % Remove control times too close to sampling times
    tl       = union(tnew, ts{experimentNumber});
    timeStep = diff(tl);
    idx      = find(timeStep < 1.0e-6);
    [~,idxtnew] = find(ismember(tnew, tl(idx)));
    tnew(idxtnew) = [];

    % Ensure that the last point coincides
    if tnew(end) ~= ts{experimentNumber}(end)
        tnew(end+1) = ts{experimentNumber}(end);
    end

    % Interpolate temperature at new control times
    Tnew = interp1(tcon_inter, T_inter, tnew, 'nearest');

    % Store result
    t_storage{experimentNumber} = tnew;
    T_storage{experimentNumber} = Tnew;
end


% Four experiments are available. Three are used for parameter estimation
% whereas the remaining is left for validation purposes. The different 
% combinations of experiments are:
switch iComb
    case{1}
        switch purpose
            case{'Est'}
                vExp = [1, 2, 3];
            case{'Val'}
                vExp = 4;
        end
    case{2}
        switch purpose
            case{'Est'}
                vExp = [1, 2, 4];
            case{'Val'}
                vExp = 3;
        end
    case{3}
        switch purpose
            case{'Est'}
                vExp = [1, 3, 4];
            case{'Val'}
                vExp = 2;
        end
    case{4}
        switch purpose
            case{'Est'}
                vExp = [2, 3, 4];
            case{'Val'}
                vExp = 1;
        end
end

% Arrange experimental data
kk = 0;
for iExp = vExp
    % Update index
    kk = kk + 1;

    % Some experiments contain several sampling times at which the variable
    % was not measured. These sampling times are removed

    %=====================================================================
    % Sh
    % Mesurements of Sh concentration
    rowsAllNaN       = all(isnan(Shewanella{iExp}), 2);
    tempTime         = ts{iExp}(~rowsAllNaN);  
    expData{kk}.tsSh = tempTime(:);                             % Sampling times [d], column format
    expData{kk}.Sh   = log10(Shewanella{iExp}(~rowsAllNaN,:));  % Measurements of Sh log10([CFU/g-fish])
    %=====================================================================

    %=====================================================================
    % Ph
    % Mesurements of Ph concentration
    rowsAllNaN       = all(isnan(Photobacterium{iExp}), 2);
    tempTime         = ts{iExp}(~rowsAllNaN);  
    expData{kk}.tsPh = tempTime(:);                                 % Sampling times [d], column format
    expData{kk}.Ph   = log10(Photobacterium{iExp}(~rowsAllNaN,:));  % Measurements of Ph log10([CFU/g-fish])
    %=====================================================================

    %=====================================================================
    % TVB
    % Mesurements of TVBN concentration
    rowsAllNaN        = all(isnan(TVBN{iExp}), 2);
    tempTime          = ts{iExp}(~rowsAllNaN);  
    expData{kk}.tsTVB = tempTime(:);                       % Sampling times [d], column format
    expData{kk}.TVB   = TVBN{iExp}(~rowsAllNaN,:)/100.0;   % Measurements of TVB-N [mg-TVB/g-fish]. Original units [mg-TVB/100 g-fish]
    %=====================================================================

    %=====================================================================
    % IMP
    % Mesurements of IMP concentration
    rowsAllNaN        = all(isnan(IMP{iExp}), 2);
    tempTime          = ts{iExp}(~rowsAllNaN);
    expData{kk}.tsIMP = tempTime(:);                % Sampling times [d], column format
    expData{kk}.IMP   = IMP{iExp}(~rowsAllNaN,:);   % Measurements of IMP
    %=====================================================================

    %=====================================================================
    % Ino
    % Mesurements of Ino concentration
    rowsAllNaN        = all(isnan(Ino{iExp}), 2);
    tempTime          = ts{iExp}(~rowsAllNaN);
    expData{kk}.tsIno = tempTime(:);                % Sampling times [d], column format
    expData{kk}.Ino   = Ino{iExp}(~rowsAllNaN,:);   % Measurements of Ino
    %=====================================================================

    %=====================================================================
    % Hx
    % Mesurements of Hx concentration
    rowsAllNaN       = all(isnan(Hx{iExp}), 2);
    tempTime         = ts{iExp}(~rowsAllNaN);
    expData{kk}.tsHx = tempTime(:);               % Sampling times [d], column format
    expData{kk}.Hx   = Hx{iExp}(~rowsAllNaN,:);   % Measurements of Hx
    %=====================================================================

    %=====================================================================
    % QIM
    % Mesurements of Hx concentration
    rowsAllNaN        = all(isnan(QIM{iExp}), 2);
    tempTime          = ts{iExp}(~rowsAllNaN);
    expData{kk}.tsQIM = tempTime(:);                % Sampling times [d], column format
    expData{kk}.QIM   = QIM{iExp}(~rowsAllNaN,:);   % Measurements of Hx
    %=====================================================================


    %=====================================================================
    % Filter temperature data to accelerate computations
    % Prepare time and temperature vectors
    t = t_storage{iExp}(:);                  % time in days
    T = T_storage{iExp}(:);                  % temperature in Celsius
    % Uniform grid except in the last sampling point
    t_filt = t(1:end-1);
    T_filt = T(1:end-1);
    % Uniform sampling rate [days]
    dt = median(diff(t_filt)); 
    Fs = 1.0/dt;                 % sampling frequency
    % Choose trend timescale (6h)
    tau = 6.0/24.0; % [days]
    fc  = 1.0/(2.0*pi*tau);   % cutoff frequency (cycles/day)
    Wn  = fc/(Fs/2.0);        % normalized cutoff
    % Butterworth low-pass filter (order 2 for stability)
    n = 2;
    [b,a] = butter(n, Wn, 'low');
    % Zero-phase filtering on uniform part
    T_lp = filtfilt(b, a, T_filt);
    % Reconstruct last value (linear interpolation/extrapolation)
    T_last = interp1(t_filt, T_lp, t(end), 'linear', 'extrap');
    % Final filtered series, same length as original
    T_filtered = [T_lp; T_last];
    % Reduce the number of sampling points
    ns      = 50;
    t_short = linspace(0,t(end),ns);
    T_short = interp1([t_filt;t(end)],T_filtered,t_short);
    % Controls
    expData{kk}.T   = T_short;       % Temperature measurements [C]
    expData{kk}.tT  = t_short;       % Times at which the temperature measurements were taken [d]
    % Slopes for the temperature measurements [C/d]
    expData{kk}.sT  = diff(expData{kk}.T)./diff(expData{kk}.tT);


end


end