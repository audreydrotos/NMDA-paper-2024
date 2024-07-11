function comparePhaseVehicle()

% Define folder paths
folderControl = "C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dmso\pre";
folderVehicle = "C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dmso\post";

% Get list of files in folder Control
control = dir(fullfile(folderControl, '*.mat'));

% Get list of files in folder DQP
Vehicle = dir(fullfile(folderVehicle, '*.mat'));

% set normalization
norm = 0;

% Create figure for histograms
figure('Position', [0 0 1400 1000]);

hold on; % Hold the plot to overlay polar plots

% Process files in folder Control
for q = 1:numel(control)
    % Read data from file in folder Control
    load(fullfile(folderControl, control(q).name));
    
    %% pull  out needed info from signal struct
    [spikeCounts synchrony phase] = analyzeAMfromStructNMDAproj(dynamicFields, 0, 0, [3:8]);

    % what modulation frequencies were presented
    stepsPerOctave = dynamicFields.fModSteps;
    freqs = 2 .^ [log2(dynamicFields.fModStart):1/dynamicFields.fModSteps:log2(dynamicFields.fModEnd)];
    logFreqs = log2(freqs);

    % make polar plot for control
    spikeRads = phase;
    magnitudes = synchrony;
    
    subplot(3, 4, q);
    for i = 1:length(spikeRads)
        % Plot each point individually and set its color from the colormap
        polarplot(spikeRads(i), logFreqs(i), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', '#e76f51');
        hold on
    end
    hold off; % Release the plot hold
    
    % Adjusting the theta direction to clockwise
    set(gca, 'ThetaDir', 'clockwise');
    
    % Displaying gridlines
    grid on 
    titleStr = control(q).name;
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)

    %Read data from file in folder DQP
    load(fullfile(folderVehicle, Vehicle(q).name));

    % analyze the data
    [spikeCounts synchrony phase] = analyzeAMfromStructNMDAproj(dynamicFields, 0, 0, [1:6]);

    % change var names
    spikeRads = phase;
    magnitudes = synchrony;
    logFreqs = log2(freqs);

    % place in correct location
    subplot(3, 4, q);
    
    % Create a polar plot
    hold on; % Hold the plot for adding points
    for i = 1:length(spikeRads)
        % Plot each point individually and set its color from the colormap
        polarplot(spikeRads(i), logFreqs(i), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', '#2a9d8f');
    end
    hold off; % Release the plot hold
    
    % Adjusting the theta direction to clockwise
    set(gca, 'ThetaDir', 'clockwise');

end