%% Script to loop through SS files to generate a summary table (output)

% Parameters
remove_ring_flag = 1;  % 1 is spindles in ring channels should be deleted
remove_ring_thres = 1.1;  % spindles with a power 1.1 times bigger in ring are deleted

%%

fileList = swa_getFiles(pwd, {'swaFile_SS'}); 


num_files = length(fileList);

load(fileList{1}, 'Info');

% pre-allocate the output structure
output = struct(...
    'frontal_ind',         [], ...
    'posterior_ind',        [], ...
    'cooccur_ind',         [], ...
    'wave_density',        [], ...
    'wave_density_frontal',  [], ...
    'wave_density_posterior',[], ...
    'wavelength',          [], ...
    'wavelength_dev',      [], ...    
    'amplitude',           [], ...
    'amplitude_dev',       [], ...
    'amplitude_thresh',    [], ...
    'power',               [], ...
    'power_dev',           [], ...
    'power_ratio',         [], ...
    'freqpeak_frontal',    [], ...
    'freqpeak_central',    [], ...
    'freqpeak_parietal',   [], ...
    'globality',           [], ...
    'globality_dev',       [], ...
    'delays',              [], ...
    'durations',           [], ...
    'topo_power',          [], ...
    'topo_power_ratio',    [], ...    
    'topo_density',        [], ...
    'topo_dens_frontal',   [], ...
    'topo_dens_posterior',  [], ...
    'topo_dens_cooccur',   [], ...
    'topo_frontal_power',  [], ...
    'topo_posterior_power', [], ...
    'topo_post_pow_rat',    [], ...
    'topo_co_power',       [], ...
    'topo_freq',           [], ...
    'time_start',          [], ...
    'participant_id',      []);

for n = 1 : num_files
    
    swa_progress_indicator('update', n, num_files);
    
    % load the swaFile
    load(fileList{n});
    
    output(n).participant_id = getfield(strsplit(Info.Recording.dataFile, '_'), {1});

    if remove_ring_flag
        ring_channels = config_ring_channels({Info.Electrodes.labels}');
        temp_data = [SS.Channels_Power];
        max_ring = max(temp_data(ring_channels,:), [], 1);
        max_inside = max(temp_data(~ring_channels,:), [], 1);

        isring = max_ring > remove_ring_thres*max_inside;
        SS(isring) = [];
        fprintf('%s found %d spindles (%d in ring were removed)\n', output(n).participant_id{1}, length(SS), sum(isring));
    end
    
    output(n).time_start = [SS.Ref_Start];

    % split by topography
    temp_data = [SS.Channels_Power];
    
    [fZones,pZones] = frontal_posterior_zones(Info.Electrodes);
    
    fPower = zeros(3, size(temp_data, 2));
    pPower = zeros(3, size(temp_data, 2));

    for ii=1:3
        fPower(ii, :) = mean(temp_data(fZones(ii,:), :), 1);
        pPower(ii, :) = mean(temp_data(pZones(ii,:), :), 1);
    end
    fPower = max(fPower,[],1);
    pPower = max(pPower,[],1);
    
    output(n).frontal_ind = fPower > 1.5 * pPower;
    output(n).posterior_ind = pPower > 1.5 * fPower;
    output(n).cooccur_ind = ~(output(n).frontal_ind | output(n).posterior_ind);
    
    % wave density (waves per minute)
    output(n).wave_density = length(SS)/(Info.Recording.dataDim(2)/Info.Recording.sRate/60);
    output(n).wave_density_frontal = sum(output(n).frontal_ind)/(Info.Recording.dataDim(2)/Info.Recording.sRate/60);
    output(n).wave_density_posterior = sum(output(n).posterior_ind)/(Info.Recording.dataDim(2)/Info.Recording.sRate/60);
    % mean wavelength
    temp_data = swa_wave_summary(SS, Info, 'wavelengths');
    output(n).wavelength = temp_data;
    output(n).wavelength_dev = std(temp_data);
    
    % mean amplitude
    temp_data = [SS.Ref_Peak2Peak];
    output(n).amplitude = temp_data;
    output(n).amplitude_dev = std(temp_data);
    output(n).amplitude_thresh = max(Info.Parameters.Ref_AmplitudeAbsolute);
    
    % mean power
    temp_data = [SS.Ref_PeakPower];
%     output(n).power = median(temp_data);
    output(n).power = temp_data;
    output(n).power_dev = std(temp_data);

    % top power ratio (top 10% of channels)
    temp_data = sort([SS.Channels_PowerRatio], 1, 'descend');
    output(n).power_ratio = mean(temp_data(1 : prctile(1:size(temp_data, 1), 10), :));
    
    % frequency peak
    % correct for spindle frequency peak using wavelet analysis
    ref_waves = [2, 5, 8];
    [peak_freq, wavelet_transform, freq_range] = ...
        swa_SS_wavelet_peak(Data, Info, SS, ref_waves, 0);

    output(n).freqpeak_frontal = peak_freq(:, 1)';
    output(n).freqpeak_central = peak_freq(:, 2)';
    output(n).freqpeak_parietal = peak_freq(:, 3)';
       
    % mean wave globality
    temp_data = swa_wave_summary(SS, Info, 'globality');
    output(n).globality = temp_data;
    output(n).globality_dev = std(temp_data);
       
    % spindle density
    temp_data = swa_wave_summary(SS, Info, 'topo_density')...
        ./(Info.Recording.dataDim(2)/Info.Recording.sRate/60);
    output(n).topo_density = temp_data;

    % power topography
    output(n).topo_power = [SS.Channels_Power];

    % get split spindle topography
    temp_data = [SS.Channels_Power];
    output(n).topo_frontal_power = temp_data(:, output(n).frontal_ind);
    output(n).topo_posterior_power = temp_data(:, output(n).posterior_ind);
    output(n).topo_co_power = temp_data(:, output(n).cooccur_ind);
    
    temp_data = [SS(output(n).posterior_ind).Channels_PowerRatio];
    output(n).topo_post_pow_rat = mean(temp_data, 2);
    
    % split density   
    output(n).topo_dens_frontal = swa_wave_summary(...
        SS(output(n).frontal_ind), Info, 'topo_density')...
        ./(Info.Recording.dataDim(2)/Info.Recording.sRate/60);

    output(n).topo_dens_posterior = swa_wave_summary(...
        SS(output(n).posterior_ind), Info, 'topo_density')...
        ./(Info.Recording.dataDim(2)/Info.Recording.sRate/60);
    
    output(n).topo_dens_cooccur = swa_wave_summary(...
        SS(output(n).cooccur_ind), Info, 'topo_density')...
        ./(Info.Recording.dataDim(2)/Info.Recording.sRate/60);
    
    % power ratio topography
    temp_data = [SS.Channels_PowerRatio];
    output(n).topo_power_ratio = temp_data;

    % frequency peak topography
    active_ind = [SS.Channels_Active];
    freqs = [SS.Channels_PeakFreq];
    freqs(active_ind) = nan;
    output(n).topo_freq = nanmean(freqs, 2);
    
    % delays & duration
    [output(n).delays, output(n).durations] = ...
        swa_SS_delays(Data, Info, SS(output(n).cooccur_ind), 0);
    
end

% Save output
save('output.mat', 'output', 'Info');
