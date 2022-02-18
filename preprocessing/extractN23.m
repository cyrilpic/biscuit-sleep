% Extract NREM stage 2 and 3 and apply basic pre-processing
% Further steps are needed until data is ready for analysis
%
% This script assumes that sleep staging is available.
%

%% Parameters

channel_location_file = '/home/mensen/Research/matlab_toolboxes/egi_257_corrected.xyz';
filter_low_cutoff = 0.5; % Hz
filter_high_cutoff = 40; % Hz

%% Load original MFF
EEG = mff_convert_to_EEGLAB();
EEG = eeg_checkset(EEG);

if size(EEG.data, 1) < 257
    EEG.data(257,:) = zeros(1,length(EEG.data));
    EEG.chanlocs = readlocs(channel_location_file);
    EEG.nbchan = length(EEG.chanlocs);
end

%% Isolate epochs of interest (produce _N23 SET file)

% Load event (sleep scoring) file
[e_file_name, e_file_path] = uigetfile('*.mat', 'Load event data file');

event_content = load(fullfile(e_file_path, e_file_name), 'event_data');
EEG.csc_event_data = event_content.event_data;
event_data = event_content.event_data;

% Extract N2N3 epochs and get rid of 4s (arousals/grob artifacts)
good_times_N23 = [];
for i=1:length(EEG.csc_event_data)
    if mod(length(good_times_N23), 2) == 0
        % Looking for new 2/3 event series
        if EEG.csc_event_data{i, 3} == 2 || EEG.csc_event_data{i, 3} == 3
            good_times_N23 = [good_times_N23 EEG.csc_event_data{i, 2}];
        elseif EEG.csc_event_data{i, 3} == 4 && EEG.csc_event_data{i-1, 3} == 4
            good_times_N23 = [good_times_N23 EEG.csc_event_data{i, 2}];
        end
    else
        % Looking for end of serie
        if EEG.csc_event_data{i, 3} ~= 2 && EEG.csc_event_data{i, 3} ~= 3
            good_times_N23 = [good_times_N23 EEG.csc_event_data{i, 2}-1/EEG.srate];
        end
    end
end
 
if mod(length(good_times_N23), 2) == 1
    % Uneven number of event
    good_times_N23 = [good_times_N23 EEG.xmax];
end
 
good_times_N23 = reshape(good_times_N23, 2, [])';

% for keeping the time (while correcting in eegrej)
EEG.orig_times = EEG.times;
 
EEG = pop_select(EEG, 'time', good_times_N23);
fprintf('N23 duration = %.2f\n', round(EEG.pnts/(EEG.srate*60*60),2))

% Save a new file with only N23
N23_file = [EEG.filename(1:end-4) '_N23.set'];
EEG = pop_saveset(EEG, 'filename', N23_file);

%% Manual artifact removal (optional)
% score artifacts with 4s

% Filter for display purposes only (manual artifact)
dEEG = pop_eegfiltnew(EEG, filter_low_cutoff, filter_high_cutoff, [], 0, [], 0);
dEEG = csc_eeg_plotter(dEEG);

EEG = pop_select(EEG, 'notime', reshape([dEEG.csc_event_data{2,2}],2,[])');
fprintf('N23 Length wo artif = %.2f\n', round(EEG.pnts/(EEG.srate*60*60),2))
clear dEEG

% save files without artifacts
N23_file_wo_artif = [EEG.filename(1:end-4) '_N23_no_artif.set'];
EEG = pop_saveset(EEG, 'filename', N23_no_artif_file);

%% Filter (produce _filtered SET file)

% Filter the data (would also remove the abrupt edges from cutting and
% sticking together from pop_select)
EEG = pop_eegfiltnew(EEG, filter_low_cutoff, filter_high_cutoff, [], 0, [], 0);

filtered_file = [N23_no_artif_file(1:end-4) '_filtered.set'];
EEG = pop_saveset(EEG, 'filename', filtered_file, 'check', 'on');

%% Detect bad channels (visually)

% Detect bad channels based on visual inspection and hide them in csc_eeg_plotter
EEG = csc_eeg_plotter(EEG);
bad_channels_manual = EEG.csc_hidden_channels;

% Apply automatic bad channel detection from EEGLAB
[~, bad_channels_auto] = pop_rejchanspec(EEG, 'freqlims',[20 40],'stdthresh',[-3.5 3.5],'plothist', 'off');

EEG.badchannels = sort(unique([bad_channels_manual bad_channels_auto]));

% Overwrite filtered_file (no changes to data only metadata)
EEG = pop_saveset(EEG, 'filename', filtered_file, 'check', 'on');

%% Reject bad channels, interpolate, and average reference (produce _nobadch_interpol_avg  SET file)

% Remove the bad channels
EEG.urchanlocs = EEG.chanlocs;
EEG = pop_select(EEG, 'nochannel', EEG.badchannels);
% Interpolate
EEG = eeg_interp(EEG, EEG.urchanlocs);
% Average reference
EEG = pop_reref(EEG, []);

% Save to disk.
nobadch_interpol_avg_file = [filtered_file(1:end-4) '_nobadch_interpol_avg'];
EEG = pop_saveset(EEG, 'filename', nobadch_interpol_avg_file, 'check', 'on');
