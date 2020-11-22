clear;

% Create figure of coverage and length
% based on Gertler & Karadi (2015) calibrated VAR simulation results

% MPM 2020-11-20


%% Settings

% File names
load_filename = fullfile('results', 'sim_var_calib.mat'); % Load results from this file
save_suffix = '.png'; % Suffix for saved figures

% Plot titles
title_vars = {'CPI', 'IP', 'EBP'}; % Plot titles for the three variables

% Specifications and legend
procs = [1 1;
         2 2;
         3 4;
         4 4]; % Index of specification (left) and associated CI (right)
legend_text = {'VAR dm', 'VAR-LA boot', 'LP-NA boot', 'LP-LA boot'}; % Legend for specifications

% Axis limits
ylim_cover = [0.4 1]; % y-limits for coverage prob plot
ylim_length = [0 10]; % y-limits for median length plot
xticks = [1 6:6:48];   % x-axis ticks

% Line specs
line_colors = [lines(2); 0.5 0.5 0.5; 0 0 0];
line_specs = {'-.', '--', ':', '-'};


%% Load results

load(load_filename);
numspec = length(specs); % No. of specifications


%% Create figure

numvars = length(settings.resp_vars); % No. of response variables
numhorz = length(settings.horzs); % No. of estimated impulse response horizons

the_f = figure('Units', 'normalize', 'Position', [0.1 0.1 0.8 0.8]);

for i=1:numvars

    % Coverage probability
    subplot(2,numvars,i);
    hold on;
    for j=1:numspec
        plot(1:numhorz, squeeze(results.coverage_prob(i,procs(j,1),:,procs(j,2))), line_specs{j}, 'Color', line_colors(j,:), 'LineWidth', 2);
    end
    the_xlim = xlim;
    plot(the_xlim, (1-settings.alpha)*[1 1], 'Color', 'k', 'LineStyle', ':'); % Nominal confidence level
    hold off;
    xlim(the_xlim);
    set(gca, 'XTick', xticks);
    xlabel('horizon');
    ylim(ylim_cover);
    title(['coverage: ', title_vars{i}]);

    % Median length
    subplot(2,numvars,numvars+i);
    hold on;
    for j=1:numspec
        plot(1:numhorz, squeeze(results.median_length(i,procs(j,1),:,procs(j,2))), line_specs{j}, 'Color', line_colors(j,:), 'LineWidth', 2);
    end
    hold off;
    set(gca, 'XTick', xticks);
    xlabel('horizon');
    ylim(ylim_length);
    title(['median length: ', title_vars{i}]);
    
    if i==numvars
        legend(legend_text, 'Location', 'SouthEast');
    end

end

% Save
status = mkdir('figures');
save_filename = fullfile('figures', 'var_calib'); % First part of file name for saved figures
if strcmp(save_suffix, '.eps')
    saveas(the_f,strcat(save_filename, save_suffix),'epsc');
else
    saveas(the_f,strcat(save_filename, save_suffix));
end

