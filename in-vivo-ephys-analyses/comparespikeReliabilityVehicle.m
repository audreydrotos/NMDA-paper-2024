function comparespikeReliabilityVehicle()

% Define folder paths
folderControl = "C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dmso\pre";
folderVehicle = "C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dmso\post";

% Get list of files in folder Control
control = dir(fullfile(folderControl, '*.mat'));

% Get list of files in folder vehicle
vehicle = dir(fullfile(folderVehicle, '*.mat'));

% Create figure for histograms
figure('Position', [0 0 1218 918]);

hold on; % Hold the plot to overlay histograms

%% make the Rate MTFs
% Process files in folder Control
for q = 1:numel(control)
    % Read data from file in folder Control
    load(fullfile(folderControl, control(q).name));

    % what modulation frequencies were presented
    stepsPerOctave = dynamicFields.fModSteps;
    freqs = 2 .^ [log2(dynamicFields.fModStart):1/dynamicFields.fModSteps:log2(dynamicFields.fModEnd)];
    logFreqs = log2(freqs);
    
    % analyze the data
    [stdevCounts] = analyzeAMspikeReliability(dynamicFields, 0, [3:8]);

    % Plot rMTF
    subplot(3, 4, q);
    
    % determind y axis label for rate MTF based on normalization method
    yLabelStr = 'stdev spike counts';

    % plot the data
    plot(log2(freqs), stdevCounts(:), 'Color', '#e76f51');
    hold on

    % label x axis ticks to show 3 tick labels
    xticks([log2(freqs(1)) log2(freqs(ceil(length(freqs)/2))) log2(freqs(end))]);
    xticklabels({num2str(freqs(1)) num2str(freqs(ceil(length(freqs)/2))) num2str(freqs(end))});
    % label x and y axes
    xlabel('mod. frequency (Hz)');
    ylabel(yLabelStr);

    % Read data from file in folder vehicle
    load(fullfile(folderVehicle, vehicle(q).name));
    
    % analyze the data
    [stdevCounts] = analyzeAMspikeReliability(dynamicFields, 0, [1:6]);
    
    % plot data on top
    subplot(3, 4, q);

    % plot the data
    plot(log2(freqs), stdevCounts(:), 'Color', '#2a9d8f');
    hold on

    ax = gca;
    titleStr = control(q).name;
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)

    % label x axis ticks to show 3 tick labels
    xticks([log2(freqs(1)) log2(freqs(ceil(length(freqs)/2))) log2(freqs(end))]);
    xticklabels({num2str(freqs(1)) num2str(freqs(ceil(length(freqs)/2))) num2str(freqs(end))});
    % label x and y axes
    xlabel('mod. frequency (Hz)');
    ylabel(yLabelStr);
    if q == 1
        legend('Control', 'Vehicle'); % Add legend
    else
    end
    hold off

    % put synchrony values in a list for later (could adjust later if want
    % outputs)
    % synchrony_control_list(q, :) = synchrony_control;
    % synchrony_vehicle_list(q, :) = synchrony_vehicle;

end

% Add labels and title
sgtitle('Spike Count Reliability Vehicle');


end %% for function