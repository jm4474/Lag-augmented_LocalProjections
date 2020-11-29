clear;

% Create figures of coverage and length
% based on AR(1) or VAR(p) simulation results

% MPM 2019-11-22


%% Settings

% DGP type
dgp_type = 'ar1_iid'; % Either 'ar1_iid', 'ar1_arch', 'ar1_homosk', or 'var_pX' (where X is integer)

% File names
load_filename = fullfile('results', strcat('sim_', dgp_type, '.mat')); % Load results from this file
save_suffix = '.png'; % Suffix for saved figures

% Plot titles
label_dm = 'delta method'; % Delta method
label_boot = {'Efron percentile', 'Hall percentile', 'Hall percentile-t'}; % Bootstrap interval types

% Axis limits
ylim_cover = [0.6 1]; % y-limits for coverage prob plot
yticks_length = -3:1:1; % y-ticks for median length plot (log10 scale)
yticklabels_length = {'0.001', '0.01', '0.1', '1', '10'}; % y-tick labels for median length plot

% Line specs
line_colors = [0 0 0; lines(7); 0.5 0.5 0.5];
line_specs = {'-o', '--o', ':o', '-.o', '-x', '--x', ':x', '-.x', '-s'};


%% Load results

load(load_filename);

numspec = length(specs); % No. of specifications

% Specification labels
specs_dm = cell(numspec,1);
specs_boot = cell(numspec,1);
for j=1:numspec
    the_lag_aug = '';
    the_bias_corr = '';
    the_har = '';
    the_bootstrap = [];
    the_boot_lag_aug = '';
    for k=1:2:length(specs{j}) % Cycle through settings
        switch specs{j}{k}
            case 'estimator'
                the_estimator = specs{j}{k+1};
            case 'lag_aug'
                if specs{j}{k+1}
                    the_lag_aug = ' LA';
                end
            case 'bias_corr'
                if specs{j}{k+1}
                    the_bias_corr = ' BC';
                end
            case 'bootstrap'
                the_bootstrap = specs{j}{k+1};
            case 'boot_lag_aug'
                if specs{j}{k+1}
                    the_boot_lag_aug = '-la';
                end
        end
    end
    specs_dm{j} = sprintf('%s%s%s%s', upper(the_estimator), the_lag_aug, the_bias_corr, the_har); % Labels for delta method specifications
    if exist('specs_p', 'var')
        specs_dm{j} = sprintf('%s %s%d', specs_dm{j}, 'p=', specs_p(j)); % VAR study: also show lag length used
    end
    if ~isempty(the_bootstrap)
        specs_boot{j} = sprintf('%s %s%s%s', specs_dm{j}, 'boot:', the_bootstrap, the_boot_lag_aug); % Labels for bootstrap specifications
    end
end

% Keep only unique entries for delta method specifications
[specs_dm, specs_dm_ind] = unique(specs_dm, 'stable');
specs_dm_ind = specs_dm_ind';

% Drop non-bootstrap specifications from bootstrap labels
empty_cells_boot = cellfun(@isempty,specs_boot);
range = 1:numspec;
specs_boot_ind = range(~empty_cells_boot);
specs_boot = specs_boot(~empty_cells_boot);


%% Create figures

status = mkdir('figures');
save_filename = fullfile('figures', dgp_type); % First part of file name for saved figures

numdgp = size(dgp.dgps, 2); % No. of DGPs
numhorz = length(settings.horzs); % No. of estimated impulse response horizons

for s=1:numdgp % For each DGP...

    for m=1:4 % For each type of confidence interval

        % Determine labels and specification indices
        if m==1 % Delta method
            the_label = label_dm;
            the_specs = specs_dm;
            the_specs_ind = specs_dm_ind;
        else % Bootstrap
            the_label = label_boot{m-1};
            the_specs = specs_boot;
            the_specs_ind = specs_boot_ind;
        end

        the_f = figure;

        % Coverage probability
        subplot(1,2,1);
        hold on;
        for j=the_specs_ind
            plot(1:numhorz, squeeze(results.coverage_prob(s,j,:,m)), line_specs{j}, 'Color', line_colors(j,:));
        end
        the_xlim = xlim;
        plot(the_xlim, (1-settings.alpha)*[1 1], 'Color', 'k', 'LineStyle', ':'); % Nominal confidence level
        hold off;
        xlim(the_xlim);
        xlabel('horizon');
        xticklabels(settings.horzs);
        ylim(ylim_cover);
        title('coverage probability');

        % Median length
        subplot(1,2,2);
        hold on;
        for j=the_specs_ind
            plot(1:numhorz, squeeze(log10(results.median_length(s,j,:,m))), line_specs{j}, 'Color', line_colors(j,:));
        end
        hold off;
        xlabel('horizon');
        xticklabels(settings.horzs);
        ylim([min(yticks_length) max(yticks_length)]);
        yticks(yticks_length);
        yticklabels(yticklabels_length);
        title('median length, log scale');
        legend(the_specs, 'Location', 'SouthEast');

        % Overall title
        sgtitle(sprintf('%s %4d%s %4.2f%s %s %s', 'T=', dgp.dgps(2,s), ', rho=', dgp.dgps(1,s), ':', the_label, 'interval'));
        
        % Save
        saveas(the_f,sprintf('%s%s%d%s%d%s', save_filename, '_dgp', s, '_intv', m, save_suffix));
        close(the_f);

    end

end




