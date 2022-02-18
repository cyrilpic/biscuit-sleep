function [ parcellation, demo ] = load_demo_parcel(participant_ids, varargin)
%LOAD_DEMO_PARCEL Load parcellatio and demo for selected participants

p = inputParser;
addRequired(p,'participant_ids');
addOptional(p,'parcellation_file','per_structure_LR_sep.csv');
addOptional(p,'demo_file','Demo.xls');
addParameter(p,'RunPCA',1,@islogical);
parse(p,participant_ids,varargin{:});

parcellation_file = p.Results.parcellation_file;
demo_file = p.Results.demo_file;

%
parcellation = readtable(parcellation_file);
parcellation.Properties.VariableNames{1} = 'participant_id';


% Remove unused participants
parcellation(~ismember(parcellation.participant_id, participant_ids), :) = [];

parcellation = sortrows(parcellation, 'participant_id');

% Add controls

controls = participant_ids(~ismember(participant_ids, parcellation.participant_id));


cont_parcel = cell(length(controls), size(parcellation, 2));
cont_parcel(:, 1) = controls;
cont_parcel(:, 2:end) = {0};

cont_parcel = cell2table(cont_parcel, 'VariableNames', parcellation.Properties.VariableNames);
parcellation = [parcellation; cont_parcel];

%% Add PCA

if p.Results.RunPCA
    X = parcellation(:, [121,127,128,130,131,132,135,136,137,138,165,166]);
    [coeff,score,~] = pca(table2array(X), 'Centered', false);

    pca_names = arrayfun(@(x) ['PCA' num2str(x)], 1:8, 'UniformOutput', false);
    pca_table = array2table(score(:, 1:8), 'VariableNames', pca_names);

    parcellation = [parcellation, pca_table];
end

%%
%load demo table
demo = readtable(demo_file);
demo.Properties.VariableNames{1} = 'participant_id';


% ensure that demo and parcellation are same and will remove accordingly
% the non-thalamic in demo table
demo(~ismember(demo.participant_id, participant_ids), :) = [];

demo(:, 'visit_number') = [];

demo = sortrows(demo, 'participant_id');

cont_demo = cell(length(controls), size(demo, 2));
cont_demo(:, 1) = controls;
cont_demo(:, 2:end) = {nan};

cont_demo = cell2table(cont_demo, 'VariableNames', demo.Properties.VariableNames);
cont_demo.gender = cellfun(@(x) x(1) == 'm', cont_demo.participant_id);
cont_demo.age = cellfun(@(x) str2num(x(2:3)), cont_demo.participant_id);
cont_demo.lesion_side = repmat({'C'}, height(cont_demo), 1);

demo = [demo; cont_demo];

lesion_side = zeros(height(demo),1);
lesion_side(strcmp(demo.lesion_side, 'L')) = 1;
lesion_side(strcmp(demo.lesion_side, 'R')) = 2;
lesion_side(strcmp(demo.lesion_side, 'LR')) = 3;
demo.lesion_side = categorical(lesion_side, 0:3, {'C', 'L', 'R', 'LR'}, 'Ordinal', 1);
demo.gender = categorical(demo.gender, [0, 1], {'f','m'});
end

