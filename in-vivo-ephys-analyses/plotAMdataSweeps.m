function plotAMdataSweeps(data, dynamicFields, sweepNumStr, signalStruct)
%% inputs to this function
i = 1; % mod frequency you want to analyze from the list
j = 1; % one modulation depth
buffer = 0.200; % 200ms on either side of the plotting
shift = 0.300*50000; % shift to make the code work

% create figure
figure('Position', [0, 0, 700, 300])

% get sweep speed & direction combos in order they were used for recording
combos = signalStruct.randomCombos;

% pull out sync onsets and offsets
syncOnsetTimes = dynamicFields.syncOnsetTimes;
syncOffsetTimes = dynamicFields.syncOffsetTimes;

% pull out spike locs
spikeLocs = dynamicFields.locs;

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

% number of times each freq and intensity combo repeated
repeats = signalStruct.repeats;

% get frequency and intensity combos in order they were used for recording
combos = signalStruct.randomCombos;

%% make blank vectors to hold the data
spikeCounts = zeros(numChannels, numFreqs, numDepths);
synchrony = zeros(numChannels, numFreqs, numDepths);

% what modulation frequencies do we want?
% find all the stim onset and offset times for the i depth and j
% frequency combo
onsets = syncOnsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));
offsets = syncOffsetTimes(combos(:,2) == modDepths(j) & combos(:,1) == freqs(i));

% just pull out the first one
onset_secs = onsets(1);
offset_secs = offsets(1);

% Make string with sweep name
sweep = append('sweep_',sweepNumStr);

% acquisiton sample rate
fs = data.header.AcquisitionSampleRate;

% name the data raw for easier access
raw = data.(sweep).analogScansHighPass;

% create a timebase
timebase = linspace(1,height(raw)/fs,height(raw));

% put onsets into seconds instead of points
onset1 = onset_secs*fs;
offset1 = offset_secs*fs;
buffer1 = buffer*fs;

% plot the data with 200ms on either side
subplot(2, 1, 1);
plot(timebase, raw, 'k')
xline(onset_secs);
xline(offset_secs);
hold on
ylim([-1.5 3.5])
xlim([onset_secs, offset_secs])
ylabel('voltage (mv)')
hold off

% plot(signalStruct.signal(round(onset2)-buffer2:round(onset2)+buffer2), 1);
%% Determine parameters for a single sweep
fs = 192000; % sample rate for playing sounds
tSound = (0:1/fs:signalStruct.soundDuration-1/fs)'; % time base for sound
sweepLength = signalStruct.sweepPeriod*fs;  % length of a single sweep in samples
startPt = signalStruct.baseDuration*fs;  % start of sound stimulus within a sweep in samples
endPt = startPt + length(tSound)-1;  % end of sound stimulus within a sweep in samples

%% Make the noise waveform.
% Set low and high frequency cutoffs for noise so that we don't drive 
% speakers at frequencies mice can't hear and that could damage speakers.
lowFreqCutoff = 4000; %Units of Hz.
highFreqCutoff = 64000; %Units of Hz.
rng(signalStruct.rngSettings); %Use the above generated seed for the random number generator. This is probably redundant, but included to be safe.
tmp = randn(1,length(tSound)); % start with wideband gaussian noise, then limit the bandwidth
B = fir1(5000,[lowFreqCutoff/(fs/2) highFreqCutoff/(fs/2)]); % compute coeff's ; 5000-order FIR filter
noise = conv(tmp,B,'same')'; % apply the filter
clear tmp;

%% Generate a vector of modulation frequencies
fMod = 2 .^ [log2(signalStruct.fModStart):1/signalStruct.fModSteps:log2(signalStruct.fModEnd)];

% Make the sound
sound = (noise.*(1+1*sin(2*pi*fMod(i)*tSound)));

subplot(2, 1, 2);
plot(tSound, sound, 'k');
hold on
xlabel('time (s)')
ylabel('sound amplitude')
hold off

%% plot everything on the same graph instead


end % for function
