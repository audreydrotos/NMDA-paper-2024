function batchAnalyzeAMdataPSTHs()

%% purpose: to analyze a folder of AM data and pull out some details about 
% the latency, etc.
% as a note: unlike the other files, there is nothing in this code to look
% at make sure the data is actually noise data. so just make sure the
% folder you are analyzing only has noise data (the other code can look for
% a particular feature like fm sweep speeds or fra intensities but that is
% harder for noise). 

% set parent directory
parentDirectory = "C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\DQP\analyze_half\pre";

% obtain a list of items in the directory
fileList = dir(fullfile(parentDirectory, '*.mat'));

% initialize matchingFiles array
matchingFiles = {};

% Loop through the fileList and filter for files starting with 'FMdata'
for i = 1:numel(fileList)
    fileName = fileList(i).name;
    if startsWith(fileName, 'AMdata')
        % Add the matching file name to the list
        matchingFiles{end+1} = fileName;
    end
end

% create empty matrices to hold variables
file_list = {}; 
inhib_index = [];
fsLatency = [];
onset_score = [];
offset_score = [];
mega_matrix = [];

% start counter for file list
j = 1;

% open figure
% figure('Position', [100 100 250 400])

% loop through list of items
for i = 1:length(matchingFiles)
    % get file name
    fileName = fullfile(parentDirectory, matchingFiles{i});

    % load structure from the file
    load(fileName);

    % analyze data
    [~, inhib_index, fsLatency, onset_score, offset_score] = AMpsthFromStruct(dynamicFields, 1);
    % 0 means no plots, 1 means plots
    hold on
    titleStr = matchingFiles{i};
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)
    subtitle(['Onset = ' num2str(onset_score) ', Offset = ' num2str(offset_score) ' Hz' ', Inhib. Index = ' num2str(inhib_index)])

    % also pull out the differences in spiking during the sound
    [spikes_pre, spikes_sound, spikes_silence, spikes_post] = compareFiringRates(dynamicFields);
    if spikes_sound > spikes_silence
        matrix = [spikes_pre spikes_sound spikes_silence spikes_post];
        mega_matrix = [mega_matrix; matrix];

        % plot(matrix,'Color', [0.5, 0.5, 0.5]);
        % hold on
        % xlim([0.5 4.5])
        % xticks([1 2 3 4])
        % xticklabels({'pre', 'sound', 'silence', 'post'})
        % ylabel('firing rate (Hz)')
    end

    % assign to structures
    file_list{j} = matchingFiles{i};
    fsLatency_list(j) = fsLatency;
    inhib_index_list(j) = inhib_index;
    onset_score_list(j) = onset_score;
    offset_score_list(j) = offset_score;

    j = j+1; % update counter

end

% avg_pre = mean(mega_matrix(:, 1));
% avg_sound = mean(mega_matrix(:, 2));
% avg_silence = mean(mega_matrix(:, 3));
% avg_post = mean(mega_matrix(:, 4));
% 
% avg_matrix = [avg_pre, avg_sound, avg_silence, avg_post];
% plot(avg_matrix, 'k', 'LineWidth', 3)
% yline(avg_pre, ':k', 'LineWidth', 2)
% hold off

% add all vars to new structure
amInfo = struct('fsLatency', fsLatency_list, ...
    'inhib_index_list', inhib_index_list, ...
    'onset_score', onset_score_list, ...
    'offset_score', offset_score_list);

end % for function