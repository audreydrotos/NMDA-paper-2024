function [spikesRepeats] = analyzeAMdataRepeats(dynamicFields)

% pull out the sweep count
sweepCount = dynamicFields.sweepCount;

% pull out the number of repeats
repeats = dynamicFields.repeats;

% how many sweeps per repeat
sweepsPerRep = sweepCount/repeats;

% spikes locs
spikeLocs = dynamicFields.locs;

% pull out sync onsets and offsets
syncOnsetTimes = dynamicFields.syncOnsetTimes;
syncOffsetTimes = dynamicFields.syncOffsetTimes;

% pull out spike locs
spikeLocs = dynamicFields.locs;

% create blank vector to hold repeats info
spikesRepeats = [];
spikeCounter = [];

% create a start index
startIndex = 1;

for i = 1:repeats
    % indexing
    endIndex = startIndex + sweepsPerRep-1;
    
    % pull out onsets and offsets
    onsets = syncOnsetTimes(startIndex:endIndex);
    offsets = syncOffsetTimes(startIndex:endIndex);

    % reset indexing vector
    startIndex = endIndex + 1;

    % blank vector to hold spikes
    spike_adj = [];

    for x = 1:sweepsPerRep
        ts = spikeLocs(spikeLocs > onsets(x) & spikeLocs < offsets(x));
        
        if ~isempty(ts)
            spikeCounter = [spikeCounter ts];
        end

    end
    
    % how many spikes in that block
    numSpikes = size(spikeCounter, 2);

    % add to master list
    spikesRepeats(i) = numSpikes;

    spikeCounter = [];

end



    





