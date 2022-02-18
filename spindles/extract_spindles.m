%% Script to go through cleaned (preprocessed) sleep EEGs and find spindles
% Locates files ending by 'cleaned.set' in the current folder and
% subfolders. Save resulting SS structure next to set file.

%% Parameters

force = 0;  % 1 overwrites existing files

% Aditional check for spindles
remove_ring_flag = 1;  % 1 is spindles in ring channels should be deleted
remove_ring_thres = 1.1;  % spindles with a power 1.1 times bigger in ring are deleted

% Filter settings
low_pass = [20, 25];
high_pass = [1, 2];

% spindle detection
detection_ref= 'grid'; % Reference method
detection_rel_ampl = [4, 2]; % Relative amplitude [high low]: Standard deviations from mean negativity
detection_method = 'power'; % power or wavelet

%% Start script

fprintf('Looking for set files underneath %s', pwd)

fileList = swa_getFiles(pwd, 'cleaned.set');
num_files = length(fileList);

fprintf('...%d found\n', num_files)

% If RAM sufficient, could be a parfor
for n = 1 : num_files

    % split the path and name
    [filePath, fileName, ext] = fileparts(fileList{n});
    
    p_id = fileName(1:6);
    save_name = fullfile(filePath, ['swaFile_SS_', p_id, '.mat']);

    if exist(save_name, 'file') && ~force
        % Skipping
        fprintf('Skipping %s\n', p_id)
        continue
    end

    % Load data
    [Data, Info] = swa_convertFromEEGLAB([fileName, ext], filePath);

    % get the default settings for spindle detection
    Info = swa_getInfoDefaults(Info, 'SS');
    
    % filter the data
    [Data.Raw] = csc_filter(Data.Raw, Info.Recording.sRate,...
        high_pass, low_pass, 8);

    % change detection parameters
    Info.Parameters.Ref_Method = detection_ref;
    Info.Parameters.Ref_AmplitudeRelative = detection_rel_ampl;
    Info.Parameters.Channels_Method = detection_method;
    
    % calculate the canonical / reference / prototypical / representative / model / illustrative wave
    [Data.SSRef, Info]  = swa_CalculateReference(Data.Raw, Info);
    
    % find the spindles in the reference
    [Data, Info, SS] = swa_FindSSRef(Data, Info);
    
    % find the waves in all channels
    [Data, Info, SS] = swa_FindSSChannels(Data, Info, SS);
    
    % Remove suspicious spindles detected at the edge of the HD-EEG net
    nSS_orig = length(SS);
    if remove_ring_flag
        ring_channels = config_ring_channels({Info.Electrodes.labels}');
        temp_data = [SS.Channels_Power];
        max_ring = max(temp_data(ring_channels,:), [], 1);
        max_inside = max(temp_data(~ring_channels,:), [], 1);

        isring = max_ring > remove_ring_thres*max_inside;
        SS(isring) = [];
        fprintf('%s found %d spindles (%s in ring were removed)\n', p_id, length(SS), length(SS)-nSS_orig);
    else
        fprintf('%s found %d spindles\n', p_id, nSS_orig);
    end
          
    % save the data
    swa_saveOutput(Data, Info, SS, save_name, 1, 0);
    
end