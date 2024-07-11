function [stdevCounts] = analyzeAMspikeReliability(signalStruct, plots, repeats)

%% purpose: analyze AM data from the structure we create
% norm = select data normalization method
%   0 = none
%   1 = normalize to max spike count
%   2 = spikes/sweep
%   3 = spikes/s

% 0 means no plots

%% pull  out needed info from signal struct

% get sweep speed & direction combos in order they were used for recording
combos = signalStruct.randomCombos;

% pull out sync onsets and offsets
syncOnsetTimes = signalStruct.syncOnsetTimes;
syncOffsetTimes = signalStruct.syncOffsetTimes;

% pull out spike locs
spikeLocs = signalStruct.locs;

% pull out pre time
preTime = signalStruct.baseDuration;
postTime = signalStruct.baseDuration + signalStruct.soundDuration; 

%% figure out when specified SAM params occurred
% what modulation frequencies were presented
stepsPerOctave = signalStruct.fModSteps;
freqs = 2 .^ [log2(signalStruct.fModStart):1/signalStruct.fModSteps:log2(signalStruct.fModEnd)];
numFreqs = length(freqs);

% make vector of all modulation depths to analyze
modDepths = signalStruct.modDepth;
numDepths = length(modDepths);

% how many channels?
numChannels = 1;

% get frequency and intensity combos in order they were used for recording
combos = signalStruct.randomCombos;

%% make blank vectors to hold the data
spikeCounts = zeros(numChannels, numFreqs, numDepths);
synchrony = zeros(numChannels, numFreqs, numDepths);
phase = zeros(numChannels, numFreqs, numDepths);

% loop through mod freq and depth combos to find number & location of
% spikes for each depth-frequency combo
for c = 1:numChannels
    for j = 1:numDepths
        for i = 1:numFreqs
            for k = 1:size(repeats, 2)
                % find all the stim onset and offset times for the i depth and j
                % frequency combo
                onsets = syncOnsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));
                offsets = syncOffsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));

                % create blank vector for spike_adj
                spike_adj = [];
                
                 % need to get spike locs by channel out of this cell array
                 spikeLocsChan = spikeLocs;
    
                 % adjust spike times relative to sync pulse onsets (stim
                 % onsets)
                 ts = spikeLocsChan(spikeLocsChan > (onsets(repeats(k))-preTime) & spikeLocsChan < (onsets(repeats(k))+postTime))-onsets(repeats(k))+preTime;
                 if ~isempty(ts > 0)
                    spike_adj = [spike_adj ts];
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
                % sweepCount = length(onsets);
                % stimDur = mean(offsets-onsets);

                numCountsReps(k) = count;

                % spikeCounts(c,i,j) = count;
                % synchrony(c,i,j) = r;
                % phase(c,i,j) = mean(spikeRads);
            end
            
            % take the standard deviation of the number of spikes at each
            % fMod
            stdevCounts(i, j) = std(numCountsReps);

        end
    end
    
end

%% make figures

% determind y axis label for rate MTF based on normalization method
yLabelStr = 'total spike counts';

if plots == 1
    % plot the figure
    figure();
    for c = 1:numChannels
        plot(log2(freqs), stdevCounts(:, :));
        hold on
        % calculate moving average
        % movMean = movmean(spikeCounts(c,:,:), stepsPerOctave);
        % plot(log2(freqs), movMean, 'LineWidth', 2, 'Color', 'k');
        hold off
        % set axis stuff
        ax = gca;
        titleStr = sprintf('Chan %d', 1+c-1);
        ax.Title.String = titleStr;
        % set y lim
        % ylim([0 1]);
        % label x axis ticks to show 3 tick labels
        xticks([log2(freqs(1)) log2(freqs(ceil(length(freqs)/2))) log2(freqs(end))]);
        xticklabels({num2str(freqs(1)) num2str(freqs(ceil(length(freqs)/2))) num2str(freqs(end))});
        % label x and y axes
        xlabel('stdev spike counts over reps')
        ylabel(yLabelStr);
    end

else 
end

end % for function
