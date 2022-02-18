function [ full_table ] = sinergia_os_create_table(output, relevant_fields, varargin)
% General a table where each line represent a single sleep marker event
% Demo and parcellation data are loaded through ``load_demo_parcel''

% Options
% ParcellationThreshold: threshold to consider that region is lesioned
% CleanTable: remove columns with only 0s
% LogVars: create log-version of relevant_fields

p = inputParser;
addRequired(p,'output',@isstruct);
addRequired(p,'relevant_fields',@iscellstr);
addParameter(p,'ParcellationThreshold',0.01,@isnumeric);
addParameter(p,'CleanTable',true,@islogical);
addParameter(p,'LogVars',true,@islogical);
parse(p,output,relevant_fields,varargin{:});

p_thres = p.Results.ParcellationThreshold;

[parcellation, demo] = load_demo_parcel(unique([output.participant_id]));
demo.group_label = categorical(demo.lesion_side ~= 'C', [0, 1], {'control', 'patient'}, 'Ordinal', 1);

% Add "categorical" lesions based on threshold
thal_struct = parcellation.Properties.VariableNames;
thal_struct = thal_struct(contains(thal_struct, 'Thal_') );
for n=1:length(thal_struct)
    if isnumeric(full_table{1,thal_struct{n}})
        parcellation.(['has_' ,thal_struct{n}]) = parcellation.(thal_struct{n}) >= pthres;
    end
end

parcellation.has_Red_N_L = parcellation.Red_N_L >= p_thres;
parcellation.has_Red_N_R = parcellation.Red_N_R >= p_thres;
parcellation.has_VTA_L = parcellation.VTA_L >= p_thres;
parcellation.has_SN_pc_L = parcellation.SN_pc_L >= p_thres;

parcellation.has_IL_MD_L = parcellation.Thal_IL_L > p_thres |parcellation.Thal_MDm_L > p_thres | parcellation.Thal_MDl_L > p_thres;
parcellation.has_IL_MD_R = parcellation.Thal_IL_R > p_thres |parcellation.Thal_MDm_R > p_thres | parcellation.Thal_MDl_R > p_thres;
parcellation.hasnt_VL_VPL_L = (parcellation.Thal_IL_L < p_thres & parcellation.Thal_IL_R <p_thres & parcellation.Thal_MDm_L < p_thres & parcellation.Thal_MDm_R < p_thres & parcellation.Thal_MDl_L < p_thres & parcellation.Thal_MDl_R < p_thres) & (parcellation.Thal_VL_L > p_thres | parcellation.Thal_VPL_L > p_thres) ;
parcellation.hasnt_VL_VPL_R = (parcellation.Thal_IL_L < p_thres & parcellation.Thal_IL_R <p_thres & parcellation.Thal_MDm_L < p_thres & parcellation.Thal_MDm_R < p_thres & parcellation.Thal_MDl_L < p_thres & parcellation.Thal_MDl_R < p_thres) & (parcellation.Thal_VL_R > p_thres | parcellation.Thal_VPL_R > p_thres) ;

parcellation.Thal_group = categorical(parcellation.has_IL_MD_L + 2*parcellation.has_IL_MD_R + 4*parcellation.hasnt_VL_VPL_R, ...
                                      [0:4], {'C', 'has_L', 'has_R', 'has_LR', 'hasnt'}, 'Ordinal', 1);

% dependent measure
n_events = cellfun(@(x) {length(x)}, {output.(relevant_fields{1})});
p_id = cellfun(@(x, y) repmat(y, 1, x), n_events, {output.participant_id}, 'UniformOutput', false);
p_id = [p_id{:}]';
full_table = table(p_id, 'VariableNames', {'participant_id'});

if isfield(output, 'night_phase')
    np = cellfun(@(x, y) repmat(y, 1, x), n_events, {output.night_phase}, 'UniformOutput', false);
    np = [np{:}]';
    to_add = table(np, 'VariableNames', {'night_phase'});
    full_table = [full_table, to_add];
end

for n = 1 : length(relevant_fields)
    
    specific_output = {output.(relevant_fields{n})};
    to_add = table(cell2mat(specific_output)', 'variableNames', relevant_fields(n));

    full_table = [full_table, to_add];
    
end
 
% ____________________
% preparing predictors
% ^^^^^^^^^^^^^^^^^^^^
% spindle location factor
frontal_indices =  {output.frontal_ind};
posterior_indices = {output.posterior_ind};
no_all = cellfun(@(x) size(x, 2), frontal_indices);

spindle_type_vector = ones(sum(no_all), 1) * 3;
for n = 1 : length(frontal_indices)
   
    % find sum of previous
    cumulative_sum = sum(no_all(1:n-1));
    
    % make relevant range for data points
    range = 1 + cumulative_sum : cumulative_sum + no_all(n);
    
    spindle_type_vector(range(frontal_indices{n})) = 1;
    spindle_type_vector(range(posterior_indices{n})) = 2;

end

% add spindle type vector to table
full_table = [full_table, table(spindle_type_vector, 'variableNames', {'spindle_type'})];
full_table.spindle_type = categorical(full_table.spindle_type, [2, 1, 3], {'posterior', 'frontal', 'co-occuring'}, 'Ordinal', 1); % assign nominal label

% parcellation and demo

combined = join(demo, parcellation);
full_table = join(full_table, combined);

if p.Results.CleanTable
    % optionnaly remove columns with only 0s (structures from
    % parcellations/demo table that would be impacted/fullfilled)
    ft_columns = full_table.Properties.VariableNames;
    for c=length(ft_columns):-1:1
        if isnumeric(full_table{:, c}) && all(full_table{:, c} == 0)
            full_table(:, c) = [];
        end
    end
end

% Add log
% use log (if necessary) to make data normally distributed
if p.Results.LogVars
    for n = 1 : length(relevant_fields)
        if contains(relevant_fields{n}, {'power', 'wavelength', 'amplitude'})
            full_table.(['log_' relevant_fields{n}]) = log(full_table.(relevant_fields{n}));
        end
    end
end

end

