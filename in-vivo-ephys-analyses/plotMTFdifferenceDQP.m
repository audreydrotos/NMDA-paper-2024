function plotMTFdifferenceDQP()

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
    
    %% pull  out needed info from signal struct
    % how many repeats?
    % repeats = dynamicFields.repeats;
    repeats = 6;
    
    % get sweep speed & direction combos in order they were used for recording
    combos = dynamicFields.randomCombos;
    
    % pull out sync onsets and offsets
    syncOnsetTimes = dynamicFields.syncOnsetTimes;
    syncOffsetTimes = dynamicFields.syncOffsetTimes;
    
    % pull out spike locs
    spikeLocs = dynamicFields.locs;
    
    % pull out pre time
    preTime = dynamicFields.baseDuration;
    postTime = dynamicFields.baseDuration + dynamicFields.soundDuration; 
    
    %% figure out when specified SAM params occurred
    % what modulation frequencies were presented
    stepsPerOctave = dynamicFields.fModSteps;
    freqs = 2 .^ [log2(dynamicFields.fModStart):1/dynamicFields.fModSteps:log2(dynamicFields.fModEnd)];
    numFreqs = length(freqs);
    
    % make vector of all modulation depths to analyze
    modDepths = dynamicFields.modDepth;
    numDepths = length(modDepths);
    
    % how many channels?
    numChannels = 1;
    
    % get frequency and intensity combos in order they were used for recording
    combos = dynamicFields.randomCombos;
    
    %% make blank vectors to hold the data
    spikeCounts_ctrl = zeros(numChannels, numFreqs, numDepths);
    synchrony_ctrl = zeros(numChannels, numFreqs, numDepths);

    spikeCounts_drug = zeros(numChannels, numFreqs, numDepths);
    synchrony_drug = zeros(numChannels, numFreqs, numDepths);
    
    % loop through mod freq and depth combos to find number & location of
    % spikes for each depth-frequency combo
    for c = 1:numChannels
        for j = 1:numDepths
            for i = 1:numFreqs
                % find all the stim onset and offset times for the i depth and j
                % frequency combo
                onsets = syncOnsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));
                offsets = syncOffsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i)) + 0.020; % add 20ms
                
                % create blank vector for spike_adj
                spike_adj = [];
                
                for k = 1+2:repeats+2 % look at last 6 repeats for control condition
                     % need to get spike locs by channel out of this cell array
                     spikeLocsChan = spikeLocs;
    
                     % adjust spike times relative to sync pulse onsets (stim
                     % onsets)
                     ts = spikeLocsChan(spikeLocsChan > (onsets(k)-preTime) & spikeLocsChan < (onsets(k)+postTime))-onsets(k)+preTime;
                     if ~isempty(ts > 0)
                        spike_adj = [spike_adj ts];
                     end
                end
    
                % calculate vector strength
                % whats the moduation frequency here
                fm = freqs(i);
    
                % convert spike times to radians
                spikeRads = mod(spike_adj, 1/fm) * 2 * pi * fm;
    
                % second, calculate vector strength. p2 of GOldberg and Brown 1969
                % method
                x = cos(spikeRads);
                y = sin(spikeRads);
                r = sqrt(sum(x)^2 + sum(y)^2)/length(x);
    
                count = length(spike_adj);
                sweepCount = length(onsets);
                stimDur = mean(offsets-onsets);
    
                spikeCounts_control(c,i,j) = count;
                synchrony_control(c,i,j) = r;
            end
        end
    end

    %%% Read data from file in folder Control
    load(fullfile(folderDQP, DQP(q).name));
     %% pull  out needed info from signal struct
    
    % get sweep speed & direction combos in order they were used for recording
    combos = dynamicFields.randomCombos;
    
    % pull out sync onsets and offsets
    syncOnsetTimes = dynamicFields.syncOnsetTimes;
    syncOffsetTimes = dynamicFields.syncOffsetTimes;
    
    % pull out spike locs
    spikeLocs = dynamicFields.locs;
    
    % pull out pre time
    preTime = dynamicFields.baseDuration;
    postTime = dynamicFields.baseDuration + dynamicFields.soundDuration; 
    
    %% figure out when specified SAM params occurred
    % what modulation frequencies were presented
    stepsPerOctave = dynamicFields.fModSteps;
    freqs = 2 .^ [log2(dynamicFields.fModStart):1/dynamicFields.fModSteps:log2(dynamicFields.fModEnd)];
    numFreqs = length(freqs);
    
    % make vector of all modulation depths to analyze
    modDepths = dynamicFields.modDepth;
    numDepths = length(modDepths);
    
    % how many channels?
    numChannels = 1;
    
    % get frequency and intensity combos in order they were used for recording
    combos = dynamicFields.randomCombos;
    
    % loop through mod freq and depth combos to find number & location of
    % spikes for each depth-frequency combo
    for c = 1:numChannels
        for j = 1:numDepths
            for i = 1:numFreqs
                % find all the stim onset and offset times for the i depth and j
                % frequency combo
                onsets = syncOnsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));
                offsets = syncOffsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i)) + 0.020;
                
                % create blank vector for spike_adj
                spike_adj = [];
                
                for k = 1:repeats
                     % need to get spike locs by channel out of this cell array
                     spikeLocsChan = spikeLocs;
    
                     % adjust spike times relative to sync pulse onsets (stim
                     % onsets)
                     ts = spikeLocsChan(spikeLocsChan > (onsets(k)-preTime) & spikeLocsChan < (onsets(k)+postTime))-onsets(k)+preTime;
                     if ~isempty(ts > 0)
                        spike_adj = [spike_adj ts];
                     end
                end
    
                % calculate vector strength
                % whats the moduation frequency here
                fm = freqs(i);
    
                % convert spike times to radians
                spikeRads = mod(spike_adj, 1/fm) * 2 * pi * fm;
    
                % second, calculate vector strength. p2 of GOldberg and Brown 1969
                % method
                x = cos(spikeRads);
                y = sin(spikeRads);
                r = sqrt(sum(x)^2 + sum(y)^2)/length(x);
    
                count = length(spike_adj);
                sweepCount = length(onsets);
                stimDur = mean(offsets-onsets);
    
                spikeCounts_drug(c,i,j) = count;
                synchrony_drug(c,i,j) = r;
            end
        end

    end


    % Plot figure
    subplot(3, 4, q);
    % make rate MTF figure
    
    % determine y axis label for rate MTF based on normalization method
    yLabelStr = 'total spike counts';

    spikeCounts_diff = spikeCounts_drug-spikeCounts_control;

    maxVal = max(spikeCounts_diff);
    minVal = min(spikeCounts_diff);

    lim = max(abs(minVal), abs(maxVal));

    % plot the data
    plot(log2(freqs), spikeCounts_diff(c,:,:), 'Color', '#2a9d8f');
    hold on
    % calculate moving average
    % movMean = movmean(spikeCounts(c,:,:), stepsPerOctave);
    % plot(log2(freqs), movMean, 'LineWidth', 2, 'Color', 'k');
    % set axis stuff
    ax = gca;
    titleStr = control(q).name;
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)
    % set y lim
    % ylim([0 1]);
    % label x axis ticks to show 3 tick labels
    xticks([log2(freqs(1)) log2(freqs(ceil(length(freqs)/2))) log2(freqs(end))]);
    xticklabels({num2str(freqs(1)) num2str(freqs(ceil(length(freqs)/2))) num2str(freqs(end))});
    % label x and y axes
    x1h = xlabel('mod. frequency (Hz)', 'Units', 'normalized');
    ylabel(yLabelStr);
    yline(0)
    ylim([-lim lim])
    hold off
end

% Add labels and title
sgtitle('Rate MTF differences DQP minus Control');

%%% make the vector strength MTFS
% Create figure for tMTFs
figure('Position', [0 0 1218 918]);

% Process files in folder Control
for q = 1:numel(control)
    % Read data from file in folder Control
    load(fullfile(folderControl, control(q).name));
    
    %% pull  out needed info from signal struct
    
    % get sweep speed & direction combos in order they were used for recording
    combos = dynamicFields.randomCombos;
    
    % pull out sync onsets and offsets
    syncOnsetTimes = dynamicFields.syncOnsetTimes;
    syncOffsetTimes = dynamicFields.syncOffsetTimes;
    
    % pull out spike locs
    spikeLocs = dynamicFields.locs;
    
    % pull out pre time
    preTime = dynamicFields.baseDuration;
    postTime = dynamicFields.baseDuration + dynamicFields.soundDuration; 
    
    %% figure out when specified SAM params occurred
    % what modulation frequencies were presented
    stepsPerOctave = dynamicFields.fModSteps;
    freqs = 2 .^ [log2(dynamicFields.fModStart):1/dynamicFields.fModSteps:log2(dynamicFields.fModEnd)];
    numFreqs = length(freqs);
    
    % make vector of all modulation depths to analyze
    modDepths = dynamicFields.modDepth;
    numDepths = length(modDepths);
    
    % how many channels?
    numChannels = 1;
    
    % number of times each freq and intensity combo repeated
    repeats = 4;
    
    % get frequency and intensity combos in order they were used for recording
    combos = dynamicFields.randomCombos;
    
    %% make blank vectors to hold the data
    synchrony_control = zeros(numChannels, numFreqs, numDepths);
    synchrony_drug = zeros(numChannels, numFreqs, numDepths);
    
    % loop through mod freq and depth combos to find number & location of
    % spikes for each depth-frequency combo
    for c = 1:numChannels
        for j = 1:numDepths
            for i = 1:numFreqs
                % find all the stim onset and offset times for the i depth and j
                % frequency combo
                onsets = syncOnsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));
                offsets = syncOffsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i)) + 0.02; % add 20 ms to account for offset response
                
                % create blank vector for spike_adj
                spike_adj = [];
                
                for k = 1+2:repeats+2
                     % need to get spike locs by channel out of this cell array
                     spikeLocsChan = spikeLocs;
    
                     % adjust spike times relative to sync pulse onsets (stim
                     % onsets)
                     ts = spikeLocsChan(spikeLocsChan > (onsets(k)-preTime) & spikeLocsChan < (onsets(k)+postTime))-onsets(k)+preTime;
                     if ~isempty(ts > 0)
                        spike_adj = [spike_adj ts];
                     end
                end
    
                % calculate vector strength
                % whats the moduation frequency here
                fm = freqs(i);
    
                % convert spike times to radians
                spikeRads = mod(spike_adj, 1/fm) * 2 * pi * fm;
    
                % second, calculate vector strength. p2 of GOldberg and Brown 1969
                % method
                x = cos(spikeRads);
                y = sin(spikeRads);
                r = sqrt(sum(x)^2 + sum(y)^2)/length(x);
    
                count = length(spike_adj);
                sweepCount = length(onsets);
                stimDur = mean(offsets-onsets);
    
                spikeCounts_control(c,i,j) = count;
                synchrony_control(c,i,j) = r;
            end
        end

    end

    %%% Read data from file in folder Control
    load(fullfile(folderDQP, DQP(q).name));
     %% pull  out needed info from signal struct
    
    % get sweep speed & direction combos in order they were used for recording
    combos = dynamicFields.randomCombos;
    
    % pull out sync onsets and offsets
    syncOnsetTimes = dynamicFields.syncOnsetTimes;
    syncOffsetTimes = dynamicFields.syncOffsetTimes;
    
    % pull out spike locs
    spikeLocs = dynamicFields.locs;
    
    % pull out pre time
    preTime = dynamicFields.baseDuration;
    postTime = dynamicFields.baseDuration + dynamicFields.soundDuration; 
    
    %% figure out when specified SAM params occurred
    % what modulation frequencies were presented
    stepsPerOctave = dynamicFields.fModSteps;
    freqs = 2 .^ [log2(dynamicFields.fModStart):1/dynamicFields.fModSteps:log2(dynamicFields.fModEnd)];
    numFreqs = length(freqs);
    
    % make vector of all modulation depths to analyze
    modDepths = dynamicFields.modDepth;
    numDepths = length(modDepths);
    
    % how many channels?
    numChannels = 1;
    
    % get frequency and intensity combos in order they were used for recording
    combos = dynamicFields.randomCombos;
    
    % loop through mod freq and depth combos to find number & location of
    % spikes for each depth-frequency combo
    for c = 1:numChannels
        for j = 1:numDepths
            for i = 1:numFreqs
                % find all the stim onset and offset times for the i depth and j
                % frequency combo
                onsets = syncOnsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));
                offsets = syncOffsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));
                
                % create blank vector for spike_adj
                spike_adj = [];
                
                for k = 1:repeats
                     % need to get spike locs by channel out of this cell array
                     spikeLocsChan = spikeLocs;
    
                     % adjust spike times relative to sync pulse onsets (stim
                     % onsets)
                     ts = spikeLocsChan(spikeLocsChan > (onsets(k)-preTime) & spikeLocsChan < (onsets(k)+postTime))-onsets(k)+preTime;
                     if ~isempty(ts > 0)
                        spike_adj = [spike_adj ts];
                     end
                end
    
                % calculate vector strength
                % whats the moduation frequency here
                fm = freqs(i);
    
                % convert spike times to radians
                spikeRads = mod(spike_adj, 1/fm) * 2 * pi * fm;
    
                % second, calculate vector strength. p2 of GOldberg and Brown 1969
                % method
                x = cos(spikeRads);
                y = sin(spikeRads);
                r = sqrt(sum(x)^2 + sum(y)^2)/length(x);
    
                count = length(spike_adj);
                sweepCount = length(onsets);
                stimDur = mean(offsets-onsets);
    
                spikeCounts_drug(c,i,j) = count;
                synchrony_drug(c,i,j) = r;
            end
        end

    end

    % Plot histogram

    synchrony_diff = synchrony_drug - synchrony_control;
    subplot(3, 4, q);
    plot(log2(freqs), synchrony_diff(c,:,:), 'Color', '#2a9d8f');
    hold on
    % set axis stuff
    ax = gca;
    ylim([-0.5 0.5])
    yline(0)
    % label x axis ticks to show 3 tick labels
    xticks([log2(freqs(1)) log2(freqs(ceil(length(freqs)/2))) log2(freqs(end))]);
    xticklabels({num2str(freqs(1)) num2str(freqs(ceil(length(freqs)/2))) num2str(freqs(end))});
    % label x and y axes
    x1h = xlabel('mod. frequency (Hz)', 'Units', 'normalized');
    ylabel('vector strength');
    titleStr = control(q).name;
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)
    hold off

end

% Add labels and title
sgtitle('Temporal MTF differences DQP minus Control');


end %% for function