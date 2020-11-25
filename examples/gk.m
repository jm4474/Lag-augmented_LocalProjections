clear;
addpath('../functions/');

% Empirical example
% Based on Gertler & Karadi (AEJ:Macro 2015), DOI: 10.1257/mac.20130329

% MPM 2020-07-11


%% Settings

% Data files from Gertler & Karadi (2015)
file_var = '../data/gk/VAR_data.csv';      % VAR series
file_iv = '../data/gk/factor_data.csv';    % External instrument series

% Specification
vars = {'logcpi', 'logip', 'ff', 'ebp', 'ff4_tc'};  % Variables in VAR
p = 12;                                             % Lag length

% Reporting
resp_var = 4;                   % Index of response variable
innov = 5;                      % Index of innovation
maxhorz = 48;                   % Maximum horizon to plot
alpha = 0.1;                    % Significance level
plot_ylim = [-4 4];             % y-axis limits
plot_xtick = [1 6:6:maxhorz];   % x-axis ticks

% Bootstrap
boot_num = 500;                     % # of repetitions
poolobj = gcp;
boot_workers = poolobj.NumWorkers;  % # of parallel workers
rng(20200711, 'twister');           % Set random number seed


%% Load data

dat = innerjoin(readtable(file_var), readtable(file_iv));
Y = dat{:,vars}; % Select variables
Y = Y(~any(isnan(Y),2),:); % Remove missing
Y = detrend(Y,0); % Remove mean
horzs = 1:maxhorz;


%% Lag-augmented local projection

disp('Lag-augmented LP: bootstrapping...');
[irs_lp_la, ~, cis_dm_lp_la, cis_boot_lp_la] = ...
    ir_estim(Y, p, horzs, ...
             'estimator', 'lp', 'lag_aug', true, ...
             'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
             'bootstrap', 'var', 'boot_num', boot_num, ...
             'boot_workers', boot_workers, 'verbose', true);

figure;
subplot(2,1,1);
plot_band(horzs, irs_lp_la, ...
          cis_dm_lp_la(1,:), cis_dm_lp_la(2,:), ...
          'Lag-augmented LP, delta method interval', ...
          plot_ylim, plot_xtick);
subplot(2,1,2);
plot_band(horzs, irs_lp_la, ...
          cis_boot_lp_la(1,:,3), cis_boot_lp_la(2,:,3), ...
          'Lag-augmented LP, percentile-t interval', ...
          plot_ylim, plot_xtick);
drawnow;


%% Non-augmented local projection

disp('Non-augmented LP: bootstrapping...');
[irs_lp_na, ~, cis_dm_lp_na, cis_boot_lp_na] = ...
    ir_estim(Y, p, horzs, ...
             'estimator', 'lp', 'lag_aug', false, ...
             'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
             'har', @(Y,bw) ewc(Y,bw), ...              % HAR estimator
             'har_bw', @(T) round(0.4*T.^(2/3)), ...    % HAR bandwidth
             'har_cv', @(bw) tinv(1-alpha/2,bw), ...    % HAR CV
             'bootstrap', 'var', 'boot_num', boot_num, ...
             'boot_workers', boot_workers, 'verbose', true);

figure;
subplot(2,1,1);
plot_band(horzs, irs_lp_na, ...
          cis_dm_lp_na(1,:), cis_dm_lp_na(2,:), ...
          'Non-augmented LP, delta method interval', ...
          plot_ylim, plot_xtick);
subplot(2,1,2);
plot_band(horzs, irs_lp_na, ...
          cis_boot_lp_na(1,:,3), cis_boot_lp_na(2,:,3), ...
          'Non-augmented LP, percentile-t interval', ...
          plot_ylim, plot_xtick);
drawnow;


%% Non-augmented VAR

disp('Non-augmented VAR: bootstrapping...');
[irs_var_na, ~, cis_dm_var_na, cis_boot_var_na] = ...
    ir_estim(Y, p, horzs, ...
             'estimator', 'var', 'lag_aug', false, ...
             'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
             'bootstrap', 'var', 'boot_num', boot_num, ...
             'boot_workers', boot_workers, 'verbose', true);

figure;
subplot(2,1,1);
plot_band(horzs, irs_var_na, ...
          cis_dm_var_na(1,:), cis_dm_var_na(2,:), ...
          'Non-augmented VAR, delta method interval', ...
          plot_ylim, plot_xtick);
subplot(2,1,2);
plot_band(horzs, irs_var_na, ...
          cis_boot_var_na(1,:,3), cis_boot_var_na(2,:,3), ...
          'Non-augmented VAR, percentile-t interval', ...
          plot_ylim, plot_xtick);
drawnow;

delete(poolobj);


%% Auxiliary plot function

function plot_band(x, y, lower, upper, plot_title, plot_ylim, plot_xtick)

    % Plot IRF and error bands

    plot(x, y, '-k', 'LineWidth', 2);
    hold on;
    plot(x, [lower; upper], 'LineStyle', '--', 'Color', 'k');
    line([min(x) max(x)], [0 0], 'Color', 'k');
    hold off;
    
    xlim([min(x) max(x)]);
    ylim(plot_ylim);
    set(gca, 'XTick', plot_xtick);
    set(gca, 'FontSize', 12);
    title(plot_title, 'FontSize', 14);

end