clear;

% Create figure of coverage and length
% based on empirically calibrated VAR simulation results

% MPM 2020-11-20


%% Settings

% Overall experiment
exper = 'gk'; % Either 'gk' (Gertler & Karadi) or 'kk' (Kilian & Kim)

% File names
load_filename = fullfile('results', strcat('sim_var_calib_', exper, '.mat')); % Load results from this file
save_filename = fullfile('figures', strcat('var_calib_', exper)); % First part of file name for saved figures
save_suffix = '.png'; % Suffix for saved figures

% Plot titles for response variables
switch exper
    case 'gk'
        title_vars = {'CPI', 'Indu. product.', 'Bond premium'};
    case 'kk'
        title_vars = {'Real activity', 'CPI inflation', 'Commod. inf.'};
end

% Specifications and legend
procs = [1 2;
         2 4]; % Index of specification (left) and associated CI (right)
legend_text = {'Lag-aug. VAR', 'Lag-aug. LP'}; % Legend for specifications

% Axis limits
switch exper
    case 'gk'
        ylim_cover = [0.75 1]; % y-limits for coverage prob plot
    case 'kk'
        ylim_cover = [0.85 1];
end
xticks = [1 12:12:48];   % x-axis ticks

% Line specs
line_colors = [lines(1); 0 0 0];
line_specs = {'--', '-'};


%% Load results

load(load_filename);
numspec = size(procs,1); % No. of specifications to plot


%% Create figure

numvar = length(settings.resp_vars); % No. of response variables
numhorz = length(settings.horzs); % No. of estimated impulse response horizons

the_f = figure('Units', 'normalize', 'Position', [0.1 0.1 0.8 0.8]);

for i=1:numvar

    % Coverage probability
    subplot(2,numvar,i);
    hold on;
    for j=1:numspec
        plot(1:numhorz, squeeze(results.coverage_prob(i,procs(j,1),:,procs(j,2))), line_specs{j}, 'Color', line_colors(j,:), 'LineWidth', 2);
    end
    plot([1 numhorz], (1-settings.alpha)*[1 1], 'Color', 'k', 'LineStyle', ':'); % Nominal confidence level
    hold off;
    set(gca, 'XTick', xticks);
    xlim([1 numhorz]);
%     xlabel('horizon');
    ylim(ylim_cover);
    title(['Coverage: ', title_vars{i}]);
    set(gca, 'FontSize', 12);

    % Median length
    subplot(2,numvar,numvar+i);
    hold on;
    for j=1:numspec
        plot(1:numhorz, squeeze(results.median_length(i,procs(j,1),:,procs(j,2))), line_specs{j}, 'Color', line_colors(j,:), 'LineWidth', 2);
    end
    hold off;
    set(gca, 'XTick', xticks);
    xlim([1 numhorz]);
%     xlabel('horizon');
    title(['Median length: ', title_vars{i}]);
    set(gca, 'FontSize', 12);
    
    if i==1
        legend(legend_text, 'Location', 'NorthWest');
    end

end

% Save
status = mkdir('figures');
if strcmp(save_suffix, '.eps')
    print(the_f,strcat(save_filename, save_suffix),'-depsc','-loose');
else
    saveas(the_f,strcat(save_filename, save_suffix));
end

