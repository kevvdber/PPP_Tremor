% cSub = "sub-POMU0A0E753F20F6332D";
% cohort = "POM";
% FirstLevel_TremorRegressor_test(cohort, cSub);

%CHECK THE FOLLOWING:
%1. MODEL NAME
%2. NUISANCE REGRESSORS
%3. THRESHOLD
%4. HIGH-PASS FILTER (INF IF COSINE USED; OTHERWISE 128)

function FirstLevel_TremorRegressor_kevvdber_05JUN2024(cohort, cSub, modelName)
%DISCUSS WITH FREEK
% WHICH TREMOR REGRESSOR (AMPLITUDE/LOG/POWER) &  'lin_unconvolved'
% 'lin_convolved' 'deriv1_unconvolved' 'deriv1_convolved' --> LOG
% (convolved vs unconvolved = of spm het doet (conditie vs regressor)
% (lin deriv = dimmer switch model)

% WHICH TIMING FOR REGRESSOR --> scans vanwege script wat regresoren maakt
% WHICH TO SET HERE matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8; treshold for masking as proportion of global (80%default) (KEVIN USED 0.1)
% WHICH TO SET HERE matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST'; (I used none in my previous script)
% HOW MANY SCANS TO EXCLUDE (5, maar dat is afhankelijk van het script wat
% je runt)
% WHY USE MULTIPLE_REG TO ADD
% matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg (omdat je dan een file
% kan toevoegen)

%% Settings & preperation
%First level name (for analysisfolder)
% modelName = "TremorRegressorModel_lowmthresh"; %4 was without ICV mask, 4_1 was with ICV mask, 4_2 is testing for hpf = -Inf; instead of hpf = Inf; 4_3 with GSR; alternative nuisance set
cNODS = 5;  %NUMBER OF DUMMYY SCANS (NODS)
TremorProjectSettings();

%% Folder management
% Make sure we are on SPM12 & Have scripts available
% rmpath(genpath(fullfile(baseHome, "common", "matlab", "spm8")));
addpath(fullfile(baseHome, "common", "matlab", "spm12"));
addpath(dirScripts);

%Find subFolder and make an analysisfolder
% cSubFolder = fullfile(dirAnalysis, "POMVisit1", "MRI", "First_Level", strcat("sub-", cSub));
cSubFolder = fullfile(dirAnalysis, "TremorNetwork", "POMVisit1", "MRI", "First_Level", cSub);
cAnalysisFolder = fullfile(cSubFolder, modelName);
if ~exist(cAnalysisFolder, 'dir')
    mkdir(cAnalysisFolder)
end

%% Get cellstring of all 3D images
cScanFile = fullfile(cSubFolder, strcat(cSub, "_ses-", cohort, "Visit1_task-rest_acq-MB8_run-1_echo-1_space-MNI152NLin6Sym_desc-preproc_bold_smoothed_6mm.nii"));
if exist(cScanFile, 'file')
    cScans = spm_select('expand', cScanFile); %4D file to 3D file
    cScans = cScans(cNODS+1:end);
    NumberOfScans = length(cScans); %Store for checks
else
    error("No preprocessed scans data for current subject");
end

%% Get tremor regressor
cTremorFile_SearchPattern = fullfile(dirTremReg, strcat(cSub, "*log.mat"));
cTremorFile_dirOutput = dir(cTremorFile_SearchPattern);
cTremorFile = strcat(cTremorFile_dirOutput(end).folder, filesep, cTremorFile_dirOutput(end).name);
cTremorData = load(cTremorFile);
allReg.(string(cTremorData.names(2))) = cTremorData.R(:,2); %lin_convolved
allReg.(string(cTremorData.names(4))) = cTremorData.R(:,4); %deriv1_convolved

%% Get nussance regressors
%Find confound regressor file
cMVTquerry = fullfile(dirAnalysis, "TremorNetwork","POMVisit1", "MRI", "First_Level", cSub, strcat(cSub, "_ses-", cohort, "Visit1_task-rest_acq-MB8_run-*_desc-confounds_timeseries2.tsv"));
cMvtDir = dir(cMVTquerry); %Retrieve file

%Use latest run path for MVT
if size(cMvtDir,1)==1
    cMvtPath=strcat(cMvtDir.folder, filesep, cMvtDir.name);
elseif isempty(cMvtDir)
    error(strcat("Querry: ", cMVTquerry, " gives zero results"));
else
    warning(strcat("Querry: ", cMVTquerry, " gives multiple results, last one used"));
    cMvtPath=strcat(cMvtDir(end).folder, filesep, cMvtDir(end).name);
end

%TremorRegressorModel
% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=["csf"; ... %cerebrospinal fluid
%     "white_matter"; ...%white matter
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     allMVTreg(eq(allMVTreg, "global_signal"))]; 

%TremorRegressorModel1
% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=[
%     strcat("a_comp_cor_", string(sprintfc('%02d', 0:7))'); ... %First 8 anatomical CompCor
%     allMVTreg(contains(allMVTreg, "cosine")); ... %All cosine's
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot"))]; %

%TremorRegressorModel2
% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=[
%     strcat("a_comp_cor_", string(sprintfc('%02d', 0:7))'); ... 
%     allMVTreg(contains(allMVTreg, "cosine")); ... 
%     allMVTreg(contains(allMVTreg, "trans")); ... 
%     allMVTreg(contains(allMVTreg, "rot")); ... 
%     allMVTreg(contains(allMVTreg, "aroma2_"))];

%TremorRegressorModel2_2
mvtReg=tdfread(cMvtPath);
allMVTreg = string(fieldnames(mvtReg));
usedRegressors=[
    "framewise_displacement"; ... 
    allMVTreg(eq(allMVTreg, "global_signal")); ...
    strcat("a_comp_cor_", string(sprintfc('%02d', 0:7))'); ... 
    allMVTreg(contains(allMVTreg, "cosine")); ... 
    allMVTreg(contains(allMVTreg, "trans")); ... 
    allMVTreg(contains(allMVTreg, "rot")); ... 
    allMVTreg(contains(allMVTreg, "aroma2_")); ...
    "std_dvars"];

%TremorRegressorModel3
% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=["csf"; ... %cerebrospinal fluid
%     "white_matter"; ...%white matter
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     allMVTreg(eq(allMVTreg, "global_signal")); ...
%     allMVTreg(contains(allMVTreg, "aroma2_"))]; 

%TremorRegressorModel4
% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=["csf"; ... %cerebrospinal fluid
%     "white_matter"; ...%white matter
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     allMVTreg(contains(allMVTreg, "aroma2_"))]; 

%TremorRegressorModel5
% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=["csf"; ... %cerebrospinal fluid
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     allMVTreg(eq(allMVTreg, "global_signal")); ...
%     allMVTreg(contains(allMVTreg, "aroma2_"))]; 

%TremorRegressorModel6
% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=["csf"; ... %cerebrospinal fluid
%     "white_matter"; ...%white matter
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     allMVTreg(eq(allMVTreg, "global_signal")); ...
%     allMVTreg(contains(allMVTreg, "aroma2_"))]; 

% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=[ "std_dvars"; ... %standardized DVARS
%     "framewise_displacement"; ... %FD
%     "csf"; ... %cerebrospinal fluid
%     "white_matter"; ...%white matter
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     allMVTreg(eq(allMVTreg, "global_signal"))]; %All aroma2 noise regressors (which are labelled as non-correlating with task-regressor)

% below the set used for model_4_1 and 4_3, the latter (extra) including
% global_signal 
% mvtReg=tdfread(cMvtPath);
% allMVTreg = string(fieldnames(mvtReg));
% usedRegressors=[ "std_dvars"; ... %standardized DVARS
%     "framewise_displacement"; ... %FD
%     strcat("a_comp_cor_", string(sprintfc('%02d', 0:7))'); ... %First 8 anatomical CompCor
%     allMVTreg(contains(allMVTreg, "cosine")); ... %All cosine's
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     allMVTreg(eq(allMVTreg, "global_signal")); ... %Global signal regression, only added for model 4_3 at june 8th, 2023
%     allMVTreg(contains(allMVTreg, "aroma2_"))]; %All aroma2 noise regressors (which are labelled as non-correlating with task-regressor)

% usedRegressors=[ "std_dvars"; ... %standardized DVARS
%     "framewise_displacement"; ... %FD
%     strcat("a_comp_cor_", string(sprintfc('%02d', 0:7))'); ... %First 8 anatomical CompCor
%     allMVTreg(contains(allMVTreg, "cosine")); ... %All cosine's
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     ]; 

% usedRegressors=[ "std_dvars"; ... %standardized DVARS
%     "framewise_displacement"; ... %FD
%     allMVTreg(contains(allMVTreg, "trans")); ... %All (12) trans motion regressors
%     allMVTreg(contains(allMVTreg, "rot")); ... %All (12) rotation motion regressors
%     ];

% usedRegressors=["trans_x"; ...
%     "trans_y"; ...
%     "trans_z"; ...
%     "rot_x"; ...
%     "rot_y"; ...
%     "rot_z"]; %All aroma2 noise regressors (which are labelled as non-correlating with task-regressor)

% User-entered Movement regressors
% condN=condN-1;
cRegCounter = 0;
for cReg = usedRegressors'
    cRegCounter = cRegCounter + 1;
    cRegValues = str2double(string(mvtReg.(cReg))); %Convert to numbers
    if isnan(cRegValues(1)); cRegValues(1)=0; end %Replace initial NaN with 0
    allReg.(cReg) = cRegValues(cNODS+1:end);
end

%OPTIONAL: add unconvolved regressors
if strcmp(modelName,'TremorRegressorModel2_5_6mm')
%     idx = numel(string(fieldnames(allReg))');
    allReg.lin_unconvolved = cTremorData.R(:,1);
    allReg.deriv1_unconvolved = cTremorData.R(:,3);
%     idx = numel(matlabbatch{1}.spm.stats.fmri_spec.sess.regress);
%     matlabbatch{1}.spm.stats.fmri_spec.sess.regress((idx+1)).name = 'lin_unconvolved';
%     matlabbatch{1}.spm.stats.fmri_spec.sess.regress((idx+2)).name = 'deriv1_unconvolved';
%     matlabbatch{1}.spm.stats.fmri_spec.sess.regress((idx+1)).val = cTremorData.R(:,1);
%     matlabbatch{1}.spm.stats.fmri_spec.sess.regress((idx+2)).val = cTremorData.R(:,3);
end

%Beta
allReg.beta = ones(size(cScans,1), 1); %add 1 for each run to the end of the matrix (beta)

% Check if all regressors have the right length
for cFieldName = string(fieldnames(allReg))'
    if size(allReg.(cFieldName),1) ~= size(cScans,1)
        error(strcat("Number of Scans and ", cFieldnName, " regressor do NOT match in size"))
    end
end

%% Create Design Matrix
%Analyses Settings
matlabbatch{1} = struct;
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(cAnalysisFolder); %outputDir
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans'; % you choose scans/secs
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.735; %TR in sec
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 8; % Microtime resolution
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1; % Microtime onset
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % Model time derivatives for HRF [1 0] for temporal derivatives
matlabbatch{1}.spm.stats.fmri_spec.volt = 1; % Do not model Volterra interactions
matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; % Global normaliziation
matlabbatch{1}.spm.stats.fmri_spec.mask = {''}; % Explicit mask
matlabbatch{1}.spm.stats.fmri_spec.mask = {'/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/masks/rmask_ICV.nii,1'}; % Explicit mask
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST'; %TR<2s gebruik FAST, anders maakt het niet (discussie centrum even navragen sara van ivan tonis groep)
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.1; %Afhankelijk van ijzerophoping --> signaalverlies, checken via 2nd level --> kevvdber: set to -Inf to include all intracranial voxels

%Session settings
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(cScans); %exclude first 5 scans
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = -Inf; %-Inf if cosine used, otherwise use default 128 seconds 
% matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; %-Inf if cosine used, otherwise use default 128 seconds 

%Add regressors the matlab way
cRegCounter = 0;
for cRegName = string(fieldnames(allReg))'
    cRegCounter = cRegCounter + 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(cRegCounter).name = char(cRegName);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(cRegCounter).val = allReg.(cRegName);
end



%Save and run the matlabbatch
cSavePath = fullfile(cAnalysisFolder, "FirstLevel_MatlabBatch_Settings.mat");
save(cSavePath, 'matlabbatch');
spm_jobman('run', matlabbatch);

%% Estimate residuals
% Settings
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = cellstr(fullfile(cAnalysisFolder, "SPM.mat"));
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

%Save and run
cSavePath = fullfile(cAnalysisFolder, "FirstLevel_MatlabBatch_Estimation.mat");
save(cSavePath, 'matlabbatch');
spm_jobman('run', matlabbatch);


%% Make contrasts
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat(1) = cellstr(fullfile(cAnalysisFolder, "SPM.mat"));
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'LIN';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 0];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'DERIV';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'lin>deriv';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'deriv>lin';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;

%Save and run
cSavePath = fullfile(cAnalysisFolder, "FirstLevel_MatlabBatch_Contrasts.mat");
save(cSavePath, 'matlabbatch');
spm_jobman('run', matlabbatch);
end