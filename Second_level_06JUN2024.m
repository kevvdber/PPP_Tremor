clear all;
load('/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/clindata/Table1OFF_LEDD_MoCA.mat');
cd('/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/Second_Level');
load('powertable4.mat');
load('Clindata_27JUN2024.mat');
load('/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/RegressorMatrix/FD.mat');
load('../Subs_PD.mat');
for i = 1:size(FD,1) FD{i,3} = nanmean(FD{i,2}); FD{i,4} = nanmedian(FD{i,2}); FD{i,5} = nanmax(FD{i,2}); end
for i = 1:size(FD,1) FD{i,6} = (FD{i,3}>=0.5); end
for i = 1:size(FD,1) FD{i,7} = (FD{i,5}>=2.4); end
for i = 1:size(FD,1) FD{i,8} = nansum(FD{i,2}>2.4); end
% power.select = power.logpower_chunk<median(power.logpower_chunk);
%Paths 
addpath('/home/common/matlab/spm12')
addpath(genpath('/project/3024011.01/TremorProject/scripts/kevin'));
addpath(genpath('/home/common/matlab/fieldtrip/qsub'));
firstleveldir = '/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/First_Level/';

Sub = cellstr(spm_select('List', fullfile(firstleveldir), 'dir', '^sub-POM.*'));
model = 'TremorRegressorModel2_4_6mm'; %model = 'PPI_antitremor';
modelsave = 'TremorRegressorModel2_4_6mm_PD_only'; modelsave = 'TremorRegressorModel2_4_6mm_PD_unflipped';

%WHICH CONTRAST
cons1 = []; cons2 = []; cons3 = []; cons4 = []; 
cons = '1'; name = 'con_0001';
cons = '2'; name = 'con_0002';
% cons = '3'; name = 'con_0003';
% cons = '4'; name = 'con_0004';
%WHICH CONTRAST

% NewTable = table;
% for i = 1:numel(SideTable.Subject)
%     for k = 1:numel(Sub)
%     if SideTable.Subject{i} == Sub{k}
%         j = size(NewTable,1)+1; 
%         NewTable(j,:) = SideTable(i,:);
%     else
%     end
%     end
% end
% SideTable = NewTable; 

% cons1 = []; cons2 = []; cons3 = []; cons4 = [];

% for i = 1:numel(SideTable.Subject)
% %     if ~gt(Table_FD.Mean_FD(i),0.5)
%     if (SideTable.Side{i}(1,1) == 'R' && FOV.FOV_include_CBLM(i)==1)
%     cSub = SideTable.Subject{i};
%     contrast1 = fullfile(firstleveldir,cSub,cTremorModel,'con_0001.nii,1');
%     contrast2 = fullfile(firstleveldir,cSub,cTremorModel,'con_0002.nii,1');
%     contrast3 = fullfile(firstleveldir,cSub,cTremorModel,'con_0003.nii,1');
%     contrast4 = fullfile(firstleveldir,cSub,cTremorModel,'con_0004.nii,1');
%     j = size(cons1,1)+1;
%     cons1{j,1} = contrast1; 
%     cons2{j,1} = contrast2; 
%     cons3{j,1} = contrast3; 
%     cons4{j,1} = contrast4; 
%     elseif (SideTable.Side{i}(1,1) == 'L' && FOV.FOV_include_CBLM(i)==1)
%     cSub = SideTable.Subject{i};
%     contrast1 = fullfile(firstleveldir,cSub,cTremorModel,'invertLR_con_0001.nii,1');
%     contrast2 = fullfile(firstleveldir,cSub,cTremorModel,'invertLR_con_0002.nii,1');
%     contrast3 = fullfile(firstleveldir,cSub,cTremorModel,'invertLR_con_0003.nii,1');
%     contrast4 = fullfile(firstleveldir,cSub,cTremorModel,'invertLR_con_0004.nii,1');
%     j = size(cons1,1)+1;
%     cons1{j,1} = contrast1; 
%     cons2{j,1} = contrast2; 
%     cons3{j,1} = contrast3; 
%     cons4{j,1} = contrast4; 
%     end
% 
% %     else
% %     end
% end

% quick check
if ~eq(sum(strcmp(FD(:,1),Sub)),numel(Sub))
    warning('No aligment');
end
if ~eq(sum(strcmp(Table1OFF.Subjects,Sub)),numel(Sub))
    warning('No aligment');
end
LEDD=[];check = [];Subj=[];
load('/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/Subs_PD.mat');
for i = 1:numel(Sub)
    if ~FD{i,6} && ~FD{i,7} && any(contains(Subs(:,1),Sub{i,1})) && ~(EyeTable.exclude(find(strcmp(EyeTable.Sub,Sub{i,1})))) %&& ~exist(fullfile(firstleveldir,Sub{i,1},model,'invert'),'dir')%&& power.RTamp_acc_POMVisit1_ON(i) > 0 && power.RTamp_acc_POMVisit1_OFF(i)>0 %~strcmp(Sub{i,1},'sub-POMU6A421CB95A963ECF') && ~strcmp(Sub{i,1},'sub-POMU98768D0FB5B292BF') %&& ~isnan(Table1OFF.LEDD(i))
           
        if ~exist(fullfile(firstleveldir,Sub{i,1},model,'invert'),'dir')
        cSub = Sub{i,1};
        contrast1 = fullfile(firstleveldir,cSub,model,'con_0001.nii,1');
        contrast2 = fullfile(firstleveldir,cSub,model,'con_0002.nii,1');
        contrast3 = fullfile(firstleveldir,cSub,model,'con_0003.nii,1');
        contrast4 = fullfile(firstleveldir,cSub,model,'con_0004.nii,1');
        elseif exist(fullfile(firstleveldir,Sub{i,1},model,'invert'),'dir')
        cSub = Sub{i,1};
        contrast1 = fullfile(firstleveldir,cSub,model,'invertLR_con_0001.nii,1');
        contrast2 = fullfile(firstleveldir,cSub,model,'invertLR_con_0002.nii,1');
        contrast3 = fullfile(firstleveldir,cSub,model,'invertLR_con_0003.nii,1');
        contrast4 = fullfile(firstleveldir,cSub,model,'invertLR_con_0004.nii,1');
        end
        j = size(cons1,1)+1;
        cons1{j,1} = contrast1; 
        cons2{j,1} = contrast2; 
        cons3{j,1} = contrast3; 
        cons4{j,1} = contrast4; 
        Subj{j,1} = Sub{i,1};
        
%         k = size(LEDD,1)+1;
%         LEDD(k,1) = Table1OFF.LEDD(i);
    else
        k = size(check,1)+1;
        check{k,1} = Sub{i,1};
        check{k,2} = FD{i,3};
        check{k,3} = FD{i,5};
        if any(contains(Subs(:,1),Sub{i,1}))
            check{k,4} = 'PD';
        elseif ~any(contains(Subs(:,1),Sub{i,1}))
            check{k,4} = 'No_PD';
        end
    end
end

%%
clear matlabbatch;
if cons == '1'
matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile('/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/Second_Level/',modelsave,'con_0001')};
scans = cons1;
if ~exist(matlabbatch{1}.spm.stats.factorial_design.dir{1}) mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1}); end
elseif cons == '2'
matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile('/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/Second_Level/',modelsave,'con_0002')};
scans = cons2;
if ~exist(matlabbatch{1}.spm.stats.factorial_design.dir{1}) mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1}); end
elseif cons == '3'
matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile('/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/Second_Level/',modelsave,'con_0003')};
scans = cons3;
if ~exist(matlabbatch{1}.spm.stats.factorial_design.dir{1}) mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1}); end
elseif cons == '4'
matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile('/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/Second_Level/',modelsave,'con_0004')};
scans = cons4;
if ~exist(matlabbatch{1}.spm.stats.factorial_design.dir{1}) mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1}); end
end

%%
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
%%
% matlabbatch{1}.spm.stats.factorial_design.cov.c = LEDD;
% matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'LEDD';
% matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
% matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
% matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/project/3024006.03/TremorProject/analyses/kevin/POMVisit1/MRI/Templates/resliced_mask_ICV.nii,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/project/3024011.01/TremorProject/analyses/kevin/TremorNetwork/POMVisit1/MRI/masks/rmask_ICV.nii,1'};
% matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/project/3024006.03/TremorProject/analyses/kevin/POMVisit1/MRI/Templates/resliced_mask_ICV.nii,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = name;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

batchdir = fullfile(matlabbatch{1}.spm.stats.factorial_design.dir{1});
cd(batchdir); 
save('matlabbatch'); 
nrun = 1; 
jobfile = {fullfile(batchdir,'matlabbatch.mat')}; 
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

