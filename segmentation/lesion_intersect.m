function lesion_intersect
% lesion_intersect calculate the intersection between structures and the lesion mask
%
% The scripts requires a given folder structure to work. Inside the `main` folder
% (defined below) there needs to be the following folder:
%  - NORMALISATION (will contain the normalised brain images and maks)
%  - PATIENTS (will contain the intersection with each structure)
%  - TABLES (will contain the results)
%
% Warning this script will produce lots of files (expect a few Gb).


%% Parameters
force = false;
main = '/path/to/MyBaseFolder';
tablename = 'per_structure_LR_sep.csv';

structure_path = '/path/to/STRUCTURES_LR/'; % Folder containing the ATLAS

% List of brain structures to consider
brain_structures = {'Precentral_L';'Precentral_R';'Frontal_Sup_2_L';'Frontal_Sup_2_R';'Frontal_Mid_2_L';'Frontal_Mid_2_R';'Frontal_Inf_Oper_L';'Frontal_Inf_Oper_R';'Frontal_Inf_Tri_L';'Frontal_Inf_Tri_R';'Frontal_Inf_Orb_2_L';'Frontal_Inf_Orb_2_R';'Rolandic_Oper_L';'Rolandic_Oper_R';'Supp_Motor_Area_L';'Supp_Motor_Area_R';'Olfactory_L';'Olfactory_R';'Frontal_Sup_Medial_L';'Frontal_Sup_Medial_R';'Frontal_Med_Orb_L';'Frontal_Med_Orb_R';'Rectus_L';'Rectus_R';'OFCmed_L';'OFCmed_R';'OFCant_L';'OFCant_R';'OFCpost_L';'OFCpost_R';'OFClat_L';'OFClat_R';'Insula_L';'Insula_R';'Cingulate_Ant_L';'Cingulate_Ant_R';'Cingulate_Mid_L';'Cingulate_Mid_R';'Cingulate_Post_L';'Cingulate_Post_R';'Hippocampus_L';'Hippocampus_R';'ParaHippocampal_L';'ParaHippocampal_R';'Amygdala_L';'Amygdala_R';'Calcarine_L';'Calcarine_R';'Cuneus_L';'Cuneus_R';'Lingual_L';'Lingual_R';'Occipital_Sup_L';'Occipital_Sup_R';'Occipital_Mid_L';'Occipital_Mid_R';'Occipital_Inf_L';'Occipital_Inf_R';'Fusiform_L';'Fusiform_R';'Postcentral_L';'Postcentral_R';'Parietal_Sup_L';'Parietal_Sup_R';'Parietal_Inf_L';'Parietal_Inf_R';'SupraMarginal_L';'SupraMarginal_R';'Angular_L';'Angular_R';'Precuneus_L';'Precuneus_R';'Paracentral_Lobule_L';'Paracentral_Lobule_R';'Caudate_L';'Caudate_R';'Putamen_L';'Putamen_R';'Pallidum_L';'Pallidum_R'; 'Heschl_L';'Heschl_R';'Temporal_Sup_L';'Temporal_Sup_R';'Temporal_Pole_Sup_L';'Temporal_Pole_Sup_R';'Temporal_Mid_L';'Temporal_Mid_R';'Temporal_Pole_Mid_L';'Temporal_Pole_Mid_R';'Temporal_Inf_L';'Temporal_Inf_R';'Cerebellum_Crus1_L';'Cerebellum_Crus1_R';'Cerebellum_Crus2_L';'Cerebellum_Crus2_R';'Cerebellum_3_L';'Cerebellum_3_R';'Cerebellum_4_5_L';'Cerebellum_4_5_R';'Cerebellum_6_L';'Cerebellum_6_R';'Cerebellum_7b_L';'Cerebellum_7b_R';'Cerebellum_8_L';'Cerebellum_8_R';'Cerebellum_9_L';'Cerebellum_9_R';'Cerebellum_10_L';'Cerebellum_10_R';'Vermis_1_2';'Vermis_3';'Vermis_4_5';'Vermis_6';'Vermis_7';'Vermis_8';'Vermis_9';'Vermis_10';'Thal_AV_L';'Thal_AV_R';'Thal_LP_L';'Thal_LP_R';'Thal_VA_L';'Thal_VA_R';'Thal_VL_L';'Thal_VL_R';'Thal_VPL_L';'Thal_VPL_R';'Thal_IL_L';'Thal_IL_R';'Thal_Re_L';'Thal_Re_R';'Thal_MDm_L';'Thal_MDm_R';'Thal_MDl_L';'Thal_MDl_R';'Thal_LGN_L';'Thal_LGN_R';'Thal_MGN_L';'Thal_MGN_R';'Thal_PuA_L';'Thal_PuA_R';'Thal_PuM_L';'Thal_PuM_R';'Thal_PuL_L';'Thal_PuL_R';'Thal_PuI_L';'Thal_PuI_R';'ACC_sub_L';'ACC_sub_R';'ACC_pre_L';'ACC_pre_R';'ACC_sup_L';'ACC_sup_R';'N_Acc_L';'N_Acc_R';'VTA_L';'VTA_R';'SN_pc_L';'SN_pc_R';'SN_pr_L';'SN_pr_R';'Red_N_L';'Red_N_R';'LC_L';'LC_R';'Raphe_D';'Raphe_M'};

% List of participants to process
% Participant are expected to be numbered starting with S and 3 digits
participant_ids = {'S001', 'S002'};

%% Prepare data

for j=1:length(participant_ids)
    % Skip if exists
    if exist([main, 'NORMALISATION/' participant_ids{j} '-seg.nii'], 'file') && ~force
        continue
    end 
    % Copy file from _5nov folder
    p_files = dir([main, '../' participant_ids{j}]);
    for k=1:numel(p_files)
        if contains(p_files(k).name, '_5nov')
            break
        end
    end
    filebase = [main, '../' participant_ids{j} '/' p_files(k).name];
    label_files = swa_getFiles(filebase,'.*-label.nii');
    if isempty(label_files)
        error(['Label not found ID=' participant_ids{j}])
    end
    copyfile(label_files{1}, [main, 'NORMALISATION/' participant_ids{j} '-seg.nii'])
    b1000_files = swa_getFiles(filebase,'.*(b0|b1000t|b1000).nii');
    if isempty(b1000_files)
        error(['B1000|B0 not found ID=' participant_ids{j}])
    end
    copyfile(b1000_files{1}, [main, 'NORMALISATION/' participant_ids{j} '-b1000.nii'])
end

%% Run normalization
b1000 = spm_select('ExtFPList', [main, 'NORMALISATION'],'^S.*(b0|b1000).*nii$');
seg = spm_select('ExtFPList', [main, 'NORMALISATION'],'^S.*seg.*nii$');

for i=1:size(b1000,1)
    bws = insertAfter(seg(i, :), 'NORMALISATION/', 'bws');
    if exist(bws(1:strfind(bws, ',')-1), 'file') && ~force
        continue
    end
    spm_normalize(b1000(i, :), seg(i, :))
end

    
function spm_normalize(b1000, seg)
    spm('defaults','fmri');
    spm_jobman('initcfg');
    matlabbatch{1}.spm.tools.MRI.MRnormseg.anat = {strip(b1000)};
    matlabbatch{1}.spm.tools.MRI.MRnormseg.les = {strip(seg)};
    matlabbatch{1}.spm.tools.MRI.MRnormseg.t2 = '';
    matlabbatch{1}.spm.tools.MRI.MRnormseg.clinicaltemplate = 1;
    matlabbatch{1}.spm.tools.MRI.MRnormseg.clean = 2;
    matlabbatch{1}.spm.tools.MRI.MRnormseg.bb = [-78 -112 -50
                                                 78 76 85];
    matlabbatch{1}.spm.tools.MRI.MRnormseg.vox = [1 1 1];
    matlabbatch{1}.spm.tools.MRI.MRnormseg.ssthresh = 0.005;
    matlabbatch{1}.spm.tools.MRI.MRnormseg.DelIntermediate = 0;
    matlabbatch{1}.spm.tools.MRI.MRnormseg.Enantiomorphic = 0;
    matlabbatch{1}.spm.tools.MRI.MRnormseg.AutoSetOrigin = 1;
    
    spm_jobman('run',matlabbatch);
end
%% Move normalization

shift = length([main, 'NORMALISATION/']);

norm_files = swa_getFiles([main, 'NORMALISATION'],'bws.*nii');
for i=1:length(norm_files)
    file = norm_files{i};
    id = file(shift+4:shift+7);
    if exist([main, 'PATIENTS/' id '/' id '.nii'], 'file') && ~force
        continue
    end
    mkdir([main, 'PATIENTS/' id])
    copyfile(file, [main, 'PATIENTS/' id '/' id '.nii'])
end

%% Intersection between structures and lesion masks

patients = [main, 'PATIENTS/'];
data = dir(patients);
is_folder = [data(:).isdir];
id = {data(is_folder).name};
id(ismember(id,{'.','..'})) = [];

% INTERSECT FILES

head = 'id,lesion,';
for i=1:length(brain_structures)
    if i < length(brain_structures)
        head = [head char(brain_structures(i)) ','];
    else
        head = [head char(brain_structures(i)) '\n'];
    end
end 

if (force) || (exist(fullfile(main, 'TABLES', tablename), 'file') == 0)
    fid=fopen(fullfile(main, 'TABLES', tablename), 'w');
    fprintf(fid, head);
    fclose(fid);
    existing_ids = {};
else
    t = readtable(fullfile(main, 'TABLES', tablename));
    existing_ids = t.id;
end

% PROCESS FILES

tic

for i=1:length(id)
    structures = {};
    c_id = char(id(i));
    if any(strcmp(existing_ids, c_id)) && ~force
        continue
    end
    % Calculate lesion volume
    volt_lesion = spm_vol(fullfile(patients, c_id, [c_id '.nii']));
    v_lesion = spm_read_vols(volt_lesion);
    volume_lesion = sum(v_lesion(:))/1000;
    structures{1} = c_id;
    structures{2} = volume_lesion;
    for j=1:length(brain_structures)
        % SPM process intersection
        spm_batch(char(id(i)), char(brain_structures(j)));
        percent_structure = intersect_struct_percentage(c_id, char(brain_structures(j)));
        structures{j+2} = percent_structure;
    end
    save_text(structures);
end

toc

% Generate intersection image in SPM
function spm_batch(id, structure)
    spm('defaults','fmri');
    spm_jobman('initcfg');

    spmbatch{1}.spm.util.imcalc.input = {
                                            fullfile(patients, id, [id '.nii,1'])
                                            fullfile(structure_path, [structure '.nii,1'])
                                            };
    spmbatch{1}.spm.util.imcalc.output = structure;
    spmbatch{1}.spm.util.imcalc.outdir = {fullfile(patients, id)};
    spmbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
    spmbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    spmbatch{1}.spm.util.imcalc.options.dmtx = 0;
    spmbatch{1}.spm.util.imcalc.options.mask = 0;
    spmbatch{1}.spm.util.imcalc.options.interp = 1;
    spmbatch{1}.spm.util.imcalc.options.dtype = 4;

    spm_jobman('run',spmbatch);
end

% Calculate the volume fraction of a lesion with respect to the size of the structure
function [percent_structure] = intersect_struct_percentage(id, structure)

    tvol_lesion = spm_vol(fullfile(patients, id, [structure '.nii']));
    vt_lesion = spm_read_vols(tvol_lesion);

    tvol_structure = spm_vol(fullfile(structure_path, [structure '.nii']));
    vt_structure = spm_read_vols(tvol_structure);

    vol_lesion = sum(vt_lesion(:))/1000;
    volume_structure = sum(vt_structure(:))/1000*prod(diag(tvol_structure.mat));
    percent_structure = vol_lesion / volume_structure;
    
    if percent_structure < 0
        percent_structure = 0;
    elseif percent_structure > 0.99
        percent_structure = 1;
    end

end

% Save line
function save_text(structures)
    
    row = [cell2mat(structures(1)) ','];
    for i=2:length(structures)
        if i < length(structures)
            row = [row num2str(cell2mat(structures(i))) ','];
        else
            row = [row num2str(cell2mat(structures(i))) '\n'];
        end
    end

    fid=fopen(fullfile(main, 'TABLES', tablename), 'a');
    fprintf(fid, row);
    fclose(fid);

end

end