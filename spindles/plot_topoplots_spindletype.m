% Plot and save for each participant topoplots of
% their power, frontal power, posterior power, and
% co-occuring power side by side

% Cell array of the participants to plot
% Names must match participant_id in the summary
% ``output'' structure

participants = {'f64', 'm57'};
normalize = false; % normalize scale across all participants
log = false;

%%
p_id = cellfun(@(x) find(strcmpi([output.participant_id], x)), participants);

if normalize
    if log
        v_min = min(arrayfun(@(x) min(mean(log(output(x).topo_power, 2))), p_id));
        v_max = max(arrayfun(@(x) max(mean(log(output(x).topo_power, 2))), p_id));
    else
        v_min = min(arrayfun(@(x) min(mean(output(x).topo_power, 2)), p_id));
        v_max = max(arrayfun(@(x) max(mean(output(x).topo_power, 2)), p_id));
    end
end

handles.fig = figure('color', 'w', ...
    'position', [100, 100, 400*4, 500]);

for i=1:numel(participants)
    clf
    j = p_id(i);
    topo_power = output(j).topo_power;
    if log
        topo_power = log(topo_power);
    end
    % General topoplot
    handles.ax(1) = subplot(1, 4, 1);
    csc_Topoplot(mean(topo_power, 2), Info.Electrodes, 'axes', handles.ax(1));
    title_text = ['Power (' participants{i} ')'];
    if log
        title_text = ['Log ', title_text];
    end
    text(0, 0.47, title_text, 'HorizontalAlignment', 'center', 'FontSize', 14)
    colorbar
    % Frontal power
    handles.ax(2) = subplot(1, 4, 2);
    csc_Topoplot(mean(topo_power(:, output(j).frontal_ind), 2), Info.Electrodes, 'axes', handles.ax(2));
    text(0, 0.47, 'Frontal', 'HorizontalAlignment', 'center', 'FontSize', 14)
    colorbar
    % Posterior power
    handles.ax(3) = subplot(1, 4, 3);
    csc_Topoplot(mean(topo_power(:, output(j).posterior_ind), 2), Info.Electrodes, 'axes', handles.ax(3));
    text(0, 0.47, 'Posterior', 'HorizontalAlignment', 'center', 'FontSize', 14)
    colorbar
    % Co-occur power
    handles.ax(4) = subplot(1, 4, 4);
    csc_Topoplot(mean(topo_power(:, output(j).cooccur_ind), 2), Info.Electrodes, 'axes', handles.ax(4));
    text(0, 0.47, 'Co-occuring', 'HorizontalAlignment', 'center', 'FontSize', 14)
    colorbar
    
    % Equalise colorbar
    if normalize
        set(handles.ax, 'clim', [v_min, v_max])
    else 
        set(handles.ax, 'clim', [min(mean(topo_power, 2)), max(mean(topo_power, 2))])
    end

    if log
        basename = ['topo_log_'];
    else
        basename = ['topo_'];
    end
    saveas(handles.fig, [basename, participants{i}, '.png']);
end
