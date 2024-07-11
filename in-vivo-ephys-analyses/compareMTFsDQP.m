function compareMTFsDQP()

% Define folder paths
folderControl = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\pre';
folderDQP = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\post';

% Get list of files in folder Control
control = dir(fullfile(folderControl, '*.mat'));

% Get list of files in folder DQP
DQP = dir(fullfile(folderDQP, '*.mat'));

% set normalization
norm = 0;

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
    [spikeCounts synchrony_control phase] = analyzeAMfromStructNMDAproj(dynamicFields, 0, 0, [3:8]);

    % Plot rMTF
    subplot(3, 4, q);
    
    % determind y axis label for rate MTF based on normalization method
    yLabelStr = 'total spike counts';

    % plot the data
    plot(log2(freqs), spikeCounts(1,:,:), 'Color', '#e76f51');
    hold on

    % label x axis ticks to show 3 tick labels
    xticks([log2(freqs(1)) log2(freqs(ceil(length(freqs)/2))) log2(freqs(end))]);
    xticklabels({num2str(freqs(1)) num2str(freqs(ceil(length(freqs)/2))) num2str(freqs(end))});
    % label x and y axes
    xlabel('mod. frequency (Hz)');
    ylabel(yLabelStr);

    % Read data from file in folder Control
    load(fullfile(folderDQP, DQP(q).name));
    
    % analyze the data
    [spikeCounts synchrony_DQP phase] = analyzeAMfromStructNMDAproj(dynamicFields, 0, 0, [1:6]);
    
    % plot data on top
    subplot(3, 4, q);

    % plot the data
    plot(log2(freqs), spikeCounts(1,:,:), 'Color', '#2a9d8f');
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
    ylim([0 inf])
    ylabel(yLabelStr);
    legend('Control', 'DQP'); % Add legend
    hold off

    % put synchrony values in a list for later
    synchrony_control_list(q, :) = synchrony_control;
    synchrony_DQP_list(q, :) = synchrony_DQP;

end

% Add labels and title
sgtitle('Rate MTFs from Control and DQP');

%%% make the vector strength MTFS
% Create figure for tMTFs
figure('Position', [0 0 1218 918]);

% Process files in folder Control
for q = 1:numel(control)
    
    % plot the data
    subplot(3, 4, q);
    % make rate MTF figure
    
    % determind y axis label for rate MTF based on normalization method
    yLabelStr = 'total spike counts';
    % plot
    subplot(3, 4, q);
    plot(log2(freqs), synchrony_control_list(q, :), 'Color', '#e76f51');
    hold on

    % plot drug condition
    subplot(3, 4, q);
    plot(log2(freqs), synchrony_DQP_list(q, :), 'Color', '#2a9d8f');
    hold on

    % set axis stuff
    ax = gca;
    % set y lim
    ylim([0 1]);

    % label x axis ticks to show 3 tick labels
    xticks([log2(freqs(1)) log2(freqs(ceil(length(freqs)/2))) log2(freqs(end))]);
    xticklabels({num2str(freqs(1)) num2str(freqs(ceil(length(freqs)/2))) num2str(freqs(end))});

    % label x and y axes
    xlabel('mod. frequency (Hz)');
    ylabel('vector strength');
    titleStr = control(q).name;
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)
    legend('Control', 'DQP'); % Add legend
    hold off
end

% Add labels and title
sgtitle('Temporal MTFs from Control and DQP');

end %% for function