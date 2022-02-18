function [data_to_filter] = csc_filter(data_to_filter, sampling_rate, high_pass, low_pass, batch_size)
% TODO: everything relevant to filtering!
% Changelog (CP):
% - Filter data by batches of channels (can help when memory low)

% NOTE: SW defaults lp: 8-12 | hp: 0.3-0.8
if nargin < 3
    low_pass = [8, 12];
    high_pass = [0.3, 0.8];
end
if nargin < 5
    batch_size = 1; % Default filters all at once
end

% design lower band-pass filter
filter_design = designfilt('bandpassiir', 'DesignMethod', 'cheby2', ...
    'StopbandFrequency1', high_pass(1), 'PassbandFrequency1', high_pass(2), ...
    'StopbandFrequency2', low_pass(2), 'PassbandFrequency2', low_pass(1), ...   
    'SampleRate', sampling_rate);

% apply filter
fprintf('Applying band pass filter from %0.1f (%0.1f) to %0.1f (%0.1f) ...\n',...
    high_pass(2), high_pass(1), low_pass(1), low_pass(2));


nbchan = size(data_to_filter, 1);
[B,~,idx] = histcounts(1:257, 'BinWidth', nbchan/batch_size);


swa_progress_indicator('inilialise', 'Filter')
for i=1:length(B)
    data_to_filter(idx==i, :) = filtfilt(filter_design, data_to_filter(idx==i, :)')';
    swa_progress_indicator('update', i, length(B))
end
