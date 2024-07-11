function compareFRsDQP()

% Define folder paths
folderControl = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\pre';
folderDQP = 'C:\Users\adrotos\Dropbox (University of Michigan)\lab\analysis\NMDAR_invivo_analysis\dqp\analyze_half\post';

% Get list of files in folder Control
control = dir(fullfile(folderControl, '*.mat'));

% Get list of files in folder DQP
DQP = dir(fullfile(folderDQP, '*.mat'));

% start counter for file list
j = 1;

% loop through list of items
for q = 1:numel(control)
    % get file name
    load(fullfile(folderControl, control(q).name));

    % pull out the repeats for the pre data
    spikesRepeats = analyzeAMdataRepeats(dynamicFields);
    
    % put some lables on the plot
    subplot(3, 4, q);
    
    % create x axis
    x = [1:8];
    plot(x, spikesRepeats, 'Color', '#e76f51')

    hold on
    titleStr = control(q).name;
    titleStr = strrep(titleStr, '_', '-');
    title(titleStr)
    xlabel('repeats')
    ylabel('num spikes')

    % read in the DQP data
    load(fullfile(folderDQP, DQP(q).name));

    % pull out the repeats for the pre data
    spikesRepeats = analyzeAMdataRepeats(dynamicFields);

    % add to plot
    subplot(3, 4, q);
    
    x2 = [9:16];
    plot(x2, spikesRepeats, 'Color', '#2a9d8f')

    if q == 1 % label the plot
        legend('Control', 'DQP');
    end

    xlim([3, 14])

    hold off % for plot

end
