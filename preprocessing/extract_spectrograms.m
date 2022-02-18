%% extract sigma power and add into table

% parameters
window_length = 10; 
selected_channel = 101;

% data structures
participants = {};
sigma = [];

%% Get all files ending with cleaned.set
fileList = swa_getFiles(pwd, 'cleaned.set');
num_files = length(fileList);

fprintf('...%d found\n', num_files)

for n = 1 : num_files
    % split the path and name
    [filePath, fileName, ext] = fileparts(fileList{n});
    
    EEG = pop_loadset(fileList{n});
    p_id = strsplit(fileName, '_');
    p_id = p_id(1);
    

    [fft, freq_range, spect_time_range, psd] = spectrogram(...
        EEG.data(selected_channel, :), ... % the data
        EEG.srate * window_length, ... % the window size
        EEG.srate * window_length/2, ... % the overlap
        EEG.srate * window_length, ... % number of points (no buffer)
        EEG.srate, ... % sampling rate
        'yaxis');
    freq_ind = freq_range < 20;
    freq_range_c = freq_range(freq_ind);
 
    % Log of the power smoothed to generate clean images
    smoothed_psd = imgaussfilt(log(psd(freq_ind, :)), 1.5);

    selected_freqs = freq_range_c > 11 & freq_range_c < 15;
    sigma_power = mean(smoothed_psd(selected_freqs, :), 1);

    sigma = [sigma, sigma_power];
    participants = [participants, repmat(p_id, 1, size(sigma_power, 2))];

end
%%
sigma_power = array2table(sigma', 'VariableNames', {'sigma'});
sigma_power.participant_id = participants';
sigma_power = movevars(sigma_power, 'participant_id', 'Before', 'sigma');
%%

% load full table to link with demographics or lesion
summary = groupsummary(full_table, 'participant_id', 'min', {'Thal_group'});
isinsummary = cellfun(@(x) any(strcmp(x, summary.participant_id)), sigma_power.participant_id);
sigma_power(~isinsummary, :) = [];
%%
sigma_table = join(sigma_power, summary);

sigma_table.sigma = double(sigma_table.sigma);

%% lme
% fitlme(sigma_table, 'sigma ~ min_Thal_group + (1 | participant_id)')
