function plotAllrMTFsDQP()

% Define folder paths
folderControl = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\pre';
folderDQP = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\post';

% Get list of files in folder Control
control = dir(fullfile(folderControl, '*.mat'));

% Get list of files in folder DQP
DQP = dir(fullfile(folderDQP, '*.mat'));

hold on; % Hold the plot to overlay polar plots

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
    load(fullfile(folderDQP, DQP(q).name));

    % analyze the data
    [spikeCounts_drug synchrony_drug phase_drug] = analyzeAMfromStructNMDAproj(dynamicFields, 0, 0, [1:6]);
    
    % what's the difference between drug and control?
    spikeCounts_diff = spikeCounts_drug - spikeCounts_ctrl;

    % normalize the difference
    spikeCounts_norm_drug = spikeCounts_diff ./ mean(spikeCounts_ctrl);

    % plot the data
    plot(log2(freqs), spikeCounts_norm_drug, 'Color', '#90a955','LineWidth', 1);

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

    % pull out all the average spike counts for each cell
    spikeDiffs_avgs_drug(q) = mean(spikeCounts_diff);
    spikeDiffs_norm_drug(q) = mean(spikeCounts_norm_drug);

    % put all the vectors together so we can average across modulation fs
    % spikeDiffs_norm_vector(:, q) = spikeCounts_diff;
    spikeDiffs_norm_vector(:, q) = spikeCounts_norm_drug;
    


end

% Define folder paths
folderControl = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dmso\pre';
folderVehicle = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dmso\post';

% Get list of files in folder Control
control = dir(fullfile(folderControl, '*.mat'));

% Get list of files in folder DQP
Vehicle = dir(fullfile(folderVehicle, '*.mat'));

hold on; % Hold the plot to overlay polar plots

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
    spikeCounts_norm_vehicle = spikeCounts_diff ./ mean(spikeCounts_ctrl);

    % plot the data
    plot(log2(freqs), spikeCounts_norm_vehicle, 'Color', '#fcbf49', 'LineWidth', 1);

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

    % pull out all the average spike counts for each cell
    spikeDiffs_avgs_vehicle(q) = mean(spikeCounts_diff);
    spikeDiffs_norm_vehicle(q) = mean(spikeCounts_norm_vehicle);

    % put all the vectors together so we can average across modulation fs
    % spikeDiffs_norm_vector_veh(:, q) = spikeCounts_diff;
    spikeDiffs_norm_vector_veh(:, q) = spikeCounts_norm_vehicle;
    
    if q == numel(control)
        % add average to the plot
        plot(log2(freqs), median(spikeDiffs_norm_vector_veh, 2)', 'Color', '#f77f00', 'LineWidth', 2)
        plot(log2(freqs), median(spikeDiffs_norm_vector, 2)', 'Color', '#31572c', 'LineWidth', 2)
    else
        continue
    end

end

% plot all the spikeCounts difs as scatter plots
figure('Position', [100 100 200 420])

% add jitter to the plot
z = size(spikeDiffs_avgs_drug, 2); % Length of the vector
x = 0.8 + 0.4 * rand(1, z);

% x = ones(size(spikeDiffs_avgs_drug, 2));
y = spikeDiffs_avgs_drug;
scatter(x, y, 'MarkerFaceColor', '#90a955', 'MarkerEdgeColor', '#90a955')
hold on

% add jitter to plot
z = size(spikeDiffs_avgs_vehicle, 2); % Length of the vector
x2 = 1.8 + 0.4 * rand(1, z);

% x2 = ones(size(spikeDiffs_avgs_vehicle, 2)) .* 2;
y = spikeDiffs_avgs_vehicle;
scatter(x2, y, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')

xlim([0.5,2.5])

xticks([1, 2]); % Specify the values where you want the ticks
xticklabels({'DQP', 'Vehicle'}); % Specify the labels for the ticks

hold off

% run wilcoxon rank-sum test for unpaired data
[p_value, h] = ranksum(spikeDiffs_avgs_drug, spikeDiffs_avgs_vehicle);
% Display the p-value and test result
fprintf('Mann-Whitney U Test on differences between drug and vehicle\n');
fprintf('P-value: %.4f\n', p_value);
fprintf('h-value: %.4f\n', h);

%% plot the normalized difference in spike Counts
% plot all the spikeCounts difs as scatter plots
figure('Position', [100 100 200 420])

% add jitter to the plot
z = size(spikeDiffs_norm_drug, 2); % Length of the vector
x = 0.8 + 0.3 * rand(1, z);

% x = ones(size(spikeDiffs_norm_drug, 2));
y = spikeDiffs_norm_drug;
scatter(x, y, 'MarkerFaceColor', '#90a955', 'MarkerEdgeColor', '#90a955')
hold on
% add average
scatter(1, median(spikeDiffs_norm_drug), 50, 'MarkerFaceColor', '#31572c', 'MarkerEdgeColor', '#31572c')

% add jitter to plot
z = size(spikeDiffs_norm_vehicle, 2); % Length of the vector
x2 = 1.8 + 0.3 * rand(1, z);

% x2 = ones(size(spikeDiffs_norm_vehicle, 2)) .* 2;
y = spikeDiffs_norm_vehicle;
scatter(x2, y, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')

% add average
scatter(2, median(spikeDiffs_norm_vehicle), 50, 'MarkerFaceColor', '#f77f00', 'MarkerEdgeColor', '#f77f00')

xlim([0.5,2.5])

xticks([1, 2]); % Specify the values where you want the ticks
xticklabels({'DQP', 'Vehicle'}); % Specify the labels for the ticks

hold off

end %% for function