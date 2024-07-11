function [spike_adj, inhib_index, fsLatency, onset_score, offset_score] = AMpsthFromStruct(dynamicFields, plots)

%% purpose:analyze AM data from a data structure
% set plots = 1 if you want plots, otherwise no plots
binSize = 0.04;

% set the parameters of the PSTH
preTime = 0.3; % time before sweep sync pulse onset to analyze
postTime = 1.3; % time after sweep sync pulse onset to analyze
latencyLimit = 0.003; % time added to pre time, then check for the first spike

% Get data sample rate
fs = dynamicFields.dataHeader.AcquisitionSampleRate;

% pull out the sound duration
soundDuration = dynamicFields.soundDuration;

% Create vector to hold time stamps corrected for sweep start times.
spike_adj = [];

% Create a vector to hold first spike latencies
fsLatency = [];

% pull out the spike locs
spikeLocs = dynamicFields.locs;

% pull out sync onsets and offsets
syncOnsetTimes = dynamicFields.syncOnsetTimes;
syncOffsetTimes = dynamicFields.syncOffsetTimes;

sweepCount = size(dynamicFields.syncOnsetTimes,2);
sweepCount = sweepCount/2; % only analyze half

% Make vector holding spike times relative to sync pulse onsets.
for i = 1:sweepCount
   
    ts = spikeLocs(spikeLocs > (syncOnsetTimes(i)-preTime) & spikeLocs < (syncOnsetTimes(i)+postTime))-syncOnsetTimes(i)+preTime;
    if ~isempty(ts)
        spike_adj = [spike_adj ts];
        
        % Find first spike after stim onset
        ts1 = min(ts(ts > preTime + latencyLimit));
        
        if ~isempty(ts1)
            fsLatency = [fsLatency ts1(1)];
        end
        
    end
    
end

% calculate 'inhibition index' (based off the DSI. -1 means inhibited by
% sound, +1 means excited by sound)

num_spikes_sound = sum(spike_adj >= preTime & spike_adj <= postTime + soundDuration);
num_spikes_quiet = sum(spike_adj < preTime | spike_adj > postTime + soundDuration);

inhib_index = (num_spikes_sound - num_spikes_quiet) / (num_spikes_sound + num_spikes_quiet);

% check if onset or offset neuron
num_spikes_onset = sum(spike_adj >= 0.05 & spike_adj < 0.1);
num_spikes_onset_compare = sum(spike_adj < 0.05 | spike_adj >= 0.1 );

num_spikes_offset = sum(spike_adj >= 0.15 & spike_adj < 0.2);
num_spikes_offset_compare = sum(spike_adj < 0.15);

% create an onset/offset score
onset_score = num_spikes_onset/num_spikes_onset_compare;
offset_score = num_spikes_offset/num_spikes_offset_compare;

fsLatency = fsLatency - preTime;

if plots == 1
    % Plot histogram
    figure()
    edges = 0:binSize:(preTime+postTime);
    norm = 0;
    if norm == 1
        h = histogram(spike_adj,edges,'Normalization','probability','FaceColor','black');
    else
        h = histogram(spike_adj,edges,'Normalization','count','FaceColor','black');
    end
    title('PSTH')
    xlabel('time (s)')
    ylabel('spike count')
    subtitle(['Inhib. Index = ' num2str(inhib_index) ', FS latency = ' num2str(mean(fsLatency)) ' s' ', Onset = ' num2str(onset_score) ', Offset = ' num2str(offset_score)])

    % plot stimulus on top
    stim = edges;
    
    % N = h.Values; % Use this if want to export the histogram counts
   
    % fprintf('First spike latency = %.2f +/- %.2f (SD) ms\n',1000*mean(fsLatency),1000*std(fsLatency))
    % fprintf('A first spike was detected within %g ms in %d of %d trials\n',1000*latencyLimit,length(fsLatency),sweepCount)
end

% create mean of fsLatency
fsLatency = mean(fsLatency); % probably we should edit this to be more like the mode of the binned data--the average
% washes out the true first spike latency. for the IC we use about 5.5 ms

end % for function
