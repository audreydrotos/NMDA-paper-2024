function compareAMdrugPSTHsDQP() 

% Define folder paths
folderControl = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\pre';
folderVehicle = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\post';

% Get list of files in folder Control
control = dir(fullfile(folderControl, '*.mat'));

% Get list of files in folder Vehicle
vehicle = dir(fullfile(folderVehicle, '*.mat'));

% Create figure for histograms
figure;

hold on; % Hold the plot to overlay histograms

% Process files in folder Control
for i = 1:numel(control)
    % Read data from file in folder Control
    load(fullfile(folderControl, control(i).name));
    
    binSize = 0.05;
    
    % set the parameters of the PSTH
    preTime = 0.25; % time before sweep sync pulse onset to analyze
    postTime = 1.25; % time after sweep sync pulse onset to analyze
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
    
    % Make vector holding spike times relative to sync pulse onsets.
    for p = sweepCount*0.25:sweepCount % anaylze 6 of the 8 repeats (the last 6)
       
        ts = spikeLocs(spikeLocs > (syncOnsetTimes(p)-preTime) & spikeLocs < (syncOnsetTimes(p)+postTime))-syncOnsetTimes(p)+preTime;

        if ~isempty(ts)
            spike_adj = [spike_adj ts];
            
            % Find first spike after stim onset
            ts1 = min(ts(ts > preTime + latencyLimit));
            
            if ~isempty(ts1)
                fsLatency = [fsLatency ts1(1)];
            end
            
        end
        
    end

    % Plot histogram
    subplot(3, 4, i);
    edges = 0:binSize:(preTime+postTime);
    norm = 0;
    if norm == 1
        h = histogram(spike_adj,edges,'Normalization','probability','FaceColor','#e76f51');
    else
        h = histogram(spike_adj,edges,'Normalization','count','FaceColor','#e76f51');
    end
    titleStr = control(i).name;
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)

    xlabel('time (s)')
    ylabel('spike count')
    

    hold on

    xline(0.25)
    xline(1.25)
    % Read data from file in folder Control
    load(fullfile(folderVehicle, vehicle(i).name));
   
    % set the parameters of the PSTH
    preTime = 0.25; % time before sweep sync pulse onset to analyze
    postTime = 1.25; % time after sweep sync pulse onset to analyze
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
    
    % Make vector holding spike times relative to sync pulse onsets.
    for p = 1:sweepCount*0.75 % anaylze 6 of the 8 repeats (the first 6)
       
        ts = spikeLocs(spikeLocs > (syncOnsetTimes(p)-preTime) & spikeLocs < (syncOnsetTimes(p)+postTime))-syncOnsetTimes(p)+preTime;
        if ~isempty(ts)
            spike_adj = [spike_adj ts];
            
            % Find first spike after stim onset
            ts1 = min(ts(ts > preTime + latencyLimit));
            
            if ~isempty(ts1)
                fsLatency = [fsLatency ts1(1)];
            end
            
        end
        
    end

    % Plot histogram
    subplot(3, 4, i);
    edges = 0:binSize:(preTime+postTime);
    norm = 0;
    if norm == 1
        h = histogram(spike_adj,edges,'Normalization','probability','FaceColor','#2a9d8f');
    else
        h = histogram(spike_adj,edges,'Normalization','count','FaceColor','#2a9d8f');
    end
    titleStr = control(i).name;
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)

    xlabel('time (s)')
    ylabel('spike count')
    if i == 1
        legend('Control', 'sound on', 'sound off', 'DQP'); % Add legend
    else
        continue
    end

    hold off 
    
end


% Add labels and title
sgtitle('Histograms from Control and DQP');