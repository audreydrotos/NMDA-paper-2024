function [file_list, AMinfo] = batchAnalyzeAMdataNMDAproj(parentDirectory, repeats)

% define parent directory
% parentDirectory = "C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\post";

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
spikeCounts_list = [];
synchrony_list = [];
spikes_pre_list = [];
spikes_post_list = [];
spikes_silence_list = [];
spikes_sound_list = [];

% start counter for file list
j = 1;

% loop through list of items
for i = 1:length(matchingFiles)
    % get file name
    fileName = fullfile(parentDirectory, matchingFiles{i});

    % load structure from the file
    load(fileName);

    % check to see if the file is an FM file
    if isfield(dynamicFields, 'modDepth')
    else
        continue
    end

    % analyze the FM data from the structure
    % this code assumes the structure is called dynamicFields
    figure();
    [spikeCounts synchrony] = analyzeAMfromStructNMDAproj(dynamicFields, 0, 0, repeats);
    % 0 for norm, also 0 for plots
    close all

    % look at changes in the firing rate
    [spikes_pre, spikes_sound, spikes_silence, spikes_post] = compareFiringRates(dynamicFields);
    
    % assign to structures
    file_list{j} = matchingFiles{i};
    spikeCounts_list(j, :) = spikeCounts;
    synchrony_list(j, :) = synchrony;
    spikes_pre_list(j) = spikes_pre;
    spikes_sound_list(j) = spikes_sound;
    spikes_silence_list(j) = spikes_silence;
    spikes_post_list(j) = spikes_post;
    j = j+1; % update counter


end

% add all vars to a new structure
AMinfo = struct('spikeCounts', spikeCounts_list, ...
    'synchrony', synchrony_list, ...
    'spikes_pre', spikes_pre_list, ...
    'spikes_sound', spikes_sound_list, ...
    'spikes_silence', spikes_silence_list, ...
    'spikes_post', spikes_post_list);

end % for function