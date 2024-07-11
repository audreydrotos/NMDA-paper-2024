function plotAllrMTFsVehicle()

% Define folder paths
folderControl = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\vehicle\analyze_half\pre';
folderVehicle = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\vehicle\analyze_half\post';

% Get list of files in folder Control
control = dir(fullfile(folderControl, '*.mat'));

% Get list of files in folder DQP
Vehicle = dir(fullfile(folderVehicle, '*.mat'));

hold on; % Hold the plot to overlay polar plots

% open a new figure
figure();

% Process files in folder Control
for q = 1:numel(control)
    % Read data from file in folder Control
    load(fullfile(folderControl, control(q).name));
    
    %% pull  out needed info from signal struct
    [spikeCounts_ctrl synchrony_ctrl phase_ctrl] = analyzeAMfromStructNMDAproj(dynamicFields, 0, 0, [3:8]);

    % what modulation frequencies were presented
    stepsPerOctave = dynamicFields.fModSteps;
    freqs = 2 .^ [log2(dynamicFields.fModStart):1/dynamicFields.fModSteps:log2(dynamicFields.fModEnd)];
    logFreqs = log2(freqs);

    %Read data from file in folder DQP
    load(fullfile(folderVehicle, Vehicle(q).name));

    % analyze the data
    [spikeCounts_drug synchrony_drug phase_drug] = analyzeAMfromStructNMDAproj(dynamicFields, 0, 0, [1:6]);
    
    % what's the difference between drug and control?
    spikeCounts_diff = spikeCounts_drug - spikeCounts_ctrl;

    % normalize the difference
    spikeCounts_norm = spikeCounts_diff ./ mean(spikeCounts_ctrl);

    % plot the data
    plot(log2(freqs), spikeCounts_norm, 'Color', '#2a9d8f');

    % set axis stuff
    ax = gca;

    % label x axis ticks to show 3 tick labels
    xticks([log2(freqs(1)) log2(freqs(ceil(length(freqs)/2))) log2(freqs(end))]);
    xticklabels({num2str(freqs(1)) num2str(freqs(ceil(length(freqs)/2))) num2str(freqs(end))});

    % label x and y axes
    xlabel('mod. frequency (Hz)');
    ylabel('spike counts');
    yline(0)
    hold on

end