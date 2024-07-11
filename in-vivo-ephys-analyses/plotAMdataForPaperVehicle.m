function plotAMdataForPaperVehicle(pre_AMinfo, post_AMinfo)

%% inputs should be pre_AMinfo and post_AMinfo
control_spikeCounts_matrix = sum(pre_AMinfo.spikeCounts, 2);
drug_spikeCounts_matrix = sum(post_AMinfo.spikeCounts, 2);

%% normalized data
figure('Position', [100 100 200 320])
y = drug_spikeCounts_matrix ./ control_spikeCounts_matrix;
x = drug_spikeCounts_matrix ./ drug_spikeCounts_matrix;
x2 = repelem(2, length(drug_spikeCounts_matrix));
scatter(x, x, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
hold on
scatter(x2, y, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
xlim([0.5,2.5])
ylim([0 2])

matrix = [x'; y'];
plot(matrix, 'Color', '#fcbf49')
% color for average f77f00
y_avg = median(y);
x_avgs = [1, 2];
y_avgs = [1, y_avg];
plot(x_avgs, y_avgs, 'Color', '#f77f00', 'LineWidth', 2);
scatter(x_avgs, y_avgs, 'MarkerFaceColor', '#f77f00', 'MarkerEdgeColor', '#f77f00')
xticks([1, 2]); % Specify the values where you want the ticks
xticklabels({'Control', 'Vehicle'}); % Specify the labels for the ticks
ylabel('normalized spike counts')
hold off

%% look at dynamic range changes normalized
figure('Position', [100 100 200 320])
control_diff_matrix = range(pre_AMinfo.spikeCounts, 2);
drug_diff_matrix = range(post_AMinfo.spikeCounts, 2);

y = control_diff_matrix ./ control_diff_matrix;
x = repelem(1, length(drug_diff_matrix));
y1 = drug_diff_matrix ./ control_diff_matrix;
x1 = repelem(2, length(drug_diff_matrix));

scatter(x, y, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
hold on
scatter(x1, y1, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
xlim([0.5,2.5])
ylim([0.2 2.2])

matrix = [y'; y1'];
plot(matrix, 'Color', '#fcbf49')
% color for average f77f00
y_avg = median(y1);
x_avgs = [1, 2];
y_avgs = [1, y_avg];
plot(x_avgs, y_avgs, 'Color', '#f77f00', 'LineWidth', 2);
scatter(x_avgs, y_avgs, 'MarkerFaceColor', '#f77f00', 'MarkerEdgeColor', '#f77f00')
xticks([1, 2]); % Specify the values where you want the ticks
xticklabels({'Control', 'Vehicle'}); % Specify the labels for the ticks
ylabel('dynamic range (spike counts)')

%% run a signed rank test for spontaneous firing rate changes (sound pre)
control = pre_AMinfo.spikes_pre;
drug = post_AMinfo.spikes_pre;

% Perform signed rank test
[p, ~, stats] = signrank(control, drug);

% Calculate median difference
median_difference = median(drug - control);

% Calculate Cohen's d
mean_difference = mean(drug - control);
std_dev_difference = std(drug - control);
cohen_d = mean_difference / std_dev_difference;

% Display the p-value
disp(['p-value: ', num2str(p)]);

% Display the test statistic
disp(['Test statistic: ', num2str(stats.signedrank)]);

% Display median difference
disp(['Median difference: ', num2str(median_difference)]);

% Display Cohen's d
disp(['Cohen''s d: ', num2str(cohen_d)]);

% plot
figure('Position', [100 100 200 320])
control_diff_matrix = pre_AMinfo.spikes_pre;
drug_diff_matrix = post_AMinfo.spikes_pre;

y = control_diff_matrix ./ control_diff_matrix;
x = repelem(1, length(control_diff_matrix));
y1 = drug_diff_matrix ./ control_diff_matrix;
x1 = repelem(2, length(control_diff_matrix));

scatter(x, y, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
hold on
scatter(x1, y1, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
xlim([0.5,2.5])
ylim([0 2])

matrix = [y; y1];
plot(matrix, 'Color', '#fcbf49')
% color for average f77f00
y_avg = median(y1);
x_avgs = [1, 2];
y_avgs = [1, y_avg];
plot(x_avgs, y_avgs, 'Color', '#f77f00', 'LineWidth', 2);
scatter(x_avgs, y_avgs, 'MarkerFaceColor', '#f77f00', 'MarkerEdgeColor', '#f77f00')
xticks([1, 2]); % Specify the values where you want the ticks
xticklabels({'Control', 'Vehicle'}); % Specify the labels for the ticks
ylabel('spontaneous firing rate (Hz)')


%% run a signed rank test for spikes during the silence
control = pre_AMinfo.spikes_silence;
drug = post_AMinfo.spikes_silence;

% Perform signed rank test
[p, ~, stats] = signrank(control, drug);

% Calculate median difference
median_difference = median(drug - control);

% Calculate Cohen's d
mean_difference = mean(drug - control);
std_dev_difference = std(drug - control);
cohen_d = mean_difference / std_dev_difference;

% Display the p-value
disp(['p-value: ', num2str(p)]);

% Display the test statistic
disp(['Test statistic: ', num2str(stats.signedrank)]);

% Display median difference
disp(['Median difference: ', num2str(median_difference)]);

% Display Cohen's d
disp(['Cohen''s d: ', num2str(cohen_d)]);

% plot
figure('Position', [100 100 200 320])
control_diff_matrix = pre_AMinfo.spikes_silence;
drug_diff_matrix = post_AMinfo.spikes_silence;

y = control_diff_matrix ./ control_diff_matrix;
x = repelem(1, length(control_diff_matrix));
y1 = drug_diff_matrix ./ control_diff_matrix;
x1 = repelem(2, length(control_diff_matrix));

scatter(x, y, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
hold on
scatter(x1, y1, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
xlim([0.5,2.5])
ylim([0 2])

matrix = [y; y1];
plot(matrix, 'Color', '#fcbf49')
% color for average f77f00
y_avg = median(y1);
x_avgs = [1, 2];
y_avgs = [1, y_avg];
plot(x_avgs, y_avgs, 'Color', '#f77f00', 'LineWidth', 2);
scatter(x_avgs, y_avgs, 'MarkerFaceColor', '#f77f00', 'MarkerEdgeColor', '#f77f00')
xticks([1, 2]); % Specify the values where you want the ticks
xticklabels({'Control', 'Vehicle'}); % Specify the labels for the ticks
ylabel('spikes in silence (Hz)')

%% run a signed rank test for dynamic range changes
control = range(pre_AMinfo.spikeCounts, 2);
drug = range(post_AMinfo.spikeCounts, 2);

% Perform signed rank test
[p, ~, stats] = signrank(control, drug);

% Calculate median difference
median_difference = median(drug - control);

% Calculate Cohen's d
mean_difference = mean(drug - control);
std_dev_difference = std(drug - control);
cohen_d = mean_difference / std_dev_difference;

% Display the p-value
disp(['p-value: ', num2str(p)]);

% Display the test statistic
disp(['Test statistic: ', num2str(stats.signedrank)]);

% Display median difference
disp(['Median difference: ', num2str(median_difference)]);

% Display Cohen's d
disp(['Cohen''s d: ', num2str(cohen_d)]);

%% run a signed rank test for spike counts
control = sum(pre_AMinfo.spikeCounts, 2);
drug = sum(post_AMinfo.spikeCounts, 2);

% Perform signed rank test
[p, ~, stats] = signrank(control, drug);

% Calculate median difference
median_difference = median(drug - control);

% Calculate Cohen's d
mean_difference = mean(drug - control);
std_dev_difference = std(drug - control);
cohen_d = mean_difference / std_dev_difference;

% Display the p-value
disp(['p-value: ', num2str(p)]);

% Display the test statistic
disp(['Test statistic: ', num2str(stats.signedrank)]);

% Display median difference
disp(['Median difference: ', num2str(median_difference)]);

% Display Cohen's d
disp(['Cohen''s d: ', num2str(cohen_d)]);

%% run a signed rank test for synchrony
control = mean(pre_AMinfo.synchrony, 2);
drug = mean(post_AMinfo.synchrony, 2);

% Perform signed rank test
[p, ~, stats] = signrank(control, drug);

% Calculate median difference
median_difference = nanmedian(drug - control);

% Calculate Cohen's d
mean_difference = nanmean(drug - control);
std_dev_difference = nanstd(drug - control);
cohen_d = mean_difference / std_dev_difference;

% Display the p-value
disp(['p-value: ', num2str(p)]);

% Display the test statistic
disp(['Test statistic: ', num2str(stats.signedrank)]);

% Display median difference
disp(['Median difference: ', num2str(median_difference)]);

% Display Cohen's d
disp(['Cohen''s d: ', num2str(cohen_d)]);

control_synchrony_matrix = nanmean(pre_AMinfo.synchrony, 2);
drug_synchrony_matrix = nanmean(post_AMinfo.synchrony, 2);

% plot
figure('Position', [100 100 200 320])

y = drug_synchrony_matrix ./ control_synchrony_matrix;
x = drug_synchrony_matrix ./ drug_synchrony_matrix;
x2 = repelem(2, length(drug_synchrony_matrix));
scatter(x, x, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
hold on
scatter(x2, y, 'MarkerFaceColor', '#fcbf49', 'MarkerEdgeColor', '#fcbf49')
xlim([0.5, 2.5])
ylim([0.5, 2.5])

matrix = [x'; y'];
plot(matrix, 'Color', '#fcbf49')
% color for average f77f00
y_avg = nanmedian(y);
x_avgs = [1, 2];
y_avgs = [1, y_avg];
plot(x_avgs, y_avgs, 'Color', '#f77f00', 'LineWidth', 2);
scatter(x_avgs, y_avgs, 'MarkerFaceColor', '#f77f00', 'MarkerEdgeColor', '#f77f00')
xticks([1, 2]); % Specify the values where you want the ticks
xticklabels({'Control', 'Vehicle'}); % Specify the labels for the ticks
ylabel('vector strength')
hold off
