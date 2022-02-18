%% heat map from lesion parcellation table

% Load lesion table
sub_table = readtable('per_structure_LR_sep.csv');
sub_table = sortrows(sub_table, 'id');

structures = strrep(sub_table.Properties.VariableNames(3:end), '_', ' ');
id = sub_table.id([3, 10, 14, 4, 8, 9, 11, 5, 7, 13, 1, 2, 6, 12, 15]);
structures_per = table2array(sub_table([3, 10, 14, 4, 8, 9, 11, 5, 7, 13, 1, 2, 6, 12, 15], 3:end));

figure('color', 'w')
axes('nextplot', 'add', ...
    'xlim', [0.5, 30.5], ...
    'ylim', [0.5, 15.5], ...
    'xTick', 1 : 30, ...
    'xTickLabelRotation', 90, ...
    'yTick', 1 : 15, ...
    'xTickLabels', structures, ...
    'yTickLabels', id, ...
    'clim', [0 , 0.5]);
imagesc(structures_per);
% axis square;
colorbar;
%colormap(flipud(hot));
%% heat map for lesion overlap


% Collect all nii in Heatmaps folder
niis = swa_getFiles('Heatmaps', '.nii');

cumulative = [];

for i=1:numel(niis)
   vol = spm_vol(niis{i});
   vt = spm_read_vols(vol);
   if isempty(cumulative)
       cumulative = vt;
   else
       cumulative = cumulative + vt;
   end
end

%%
vol_hm = vol;
vol_hm.descrip = 'Heatmap';
vol_hm.fname = 'Heatmaps/heatmap.nii';

vol_hm = spm_create_vol(vol_hm);
spm_write_vol(vol_hm, cumulative);
