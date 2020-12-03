clear;

% Create table of coverage and length
% based on AR(1) and bivariate VAR(p) simulation results


%% Settings

% DGP type
dgp_type = 'ar1_iid'; % Either 'ar1_iid', 'ar1_arch', 'ar1_homosk', or 'var_pX' (where X is integer)

% File names
load_filename = fullfile('results', strcat('sim_', dgp_type, '.mat'));  % Load results from this file
save_filename = dgp_type; % File name for saved table

% Select DGPs
T_sel = 240;                            % Sample size (single number)
switch dgp_type(1:3)
    case 'ar1'
        rhos_sel = [0.0 0.5 0.95 1];    % rho values (array)
    case 'var'
        rhos_sel = [0.0 0.5 0.95 1];
end

% Select CI procedures
switch dgp_type(1:3)
    case 'ar1'
        proc_names = {'$\text{LP-LA}_b$', '$\text{LP-LA}$', '$\text{LP}_b$', '$\text{LP}$', '$\text{AR-LA}_b$', '$\text{AR}$'};
        procs = [7 4; % First index: inference procedure; second index: type of confidence interval
                 7 1;
                 6 4;
                 6 1
                 4 2
                 1 1];
    case 'var'
        proc_names = {'$\text{LP-LA}_b$', '$\text{LP-LA}_b^8$', '$\text{LP}_b$', '$\text{AR-LA}_b$', '$\text{AR}$'};
        procs = [4 4;
                 5 4;
                 3 4;
                 2 2;
                 1 1];
end


%% Load results

load(load_filename);

% Pick out indices of selected DGPs
numdgp_sel = length(rhos_sel);
dgp_sel = zeros(1,numdgp_sel);
for j=1:numdgp_sel
    dgp_sel(j) = find(dgp.dgps(1,:)==rhos_sel(j) & dgp.dgps(2,:)==T_sel);
end

numproc = length(proc_names); % Number of inference procedures


%% Write table

status = mkdir('tables');
f = fopen(fullfile('tables', strcat(save_filename, '.tex')), 'w'); % Open file for writing

fprintf(f, '%s%s%s%s%s\n', '\begin{tabular}{r|', repmat('c', 1, numproc), '|', repmat('c', 1, numproc), '}');
fprintf(f, '%s%d%s%d%s\n', '& \multicolumn{', numproc, '}{c|}{Coverage} & \multicolumn{', numproc, '}{c}{Median length} \\');
fprintf(f, '%s', '$h$');
for i=1:2
    for j=1:numproc
        fprintf(f, '%s%s', ' & ', proc_names{j});
    end
end
fprintf(f, '%s\n%s\n', ' \\', '\hline');

for d=dgp_sel
    rho = dgp.dgps(1,d);
    fprintf(f, '%s%d%s%4.2f%s\n', '\multicolumn{', 1+2*numproc, '}{c}{$\rho = ', rho, '$} \\');
    for ih=1:length(settings.horzs)
        h = settings.horzs(ih);
        fprintf(f, '%3d', h);
        for j=1:numproc
            cp = results.coverage_prob(d,procs(j,1),ih,procs(j,2));
            fprintf(f, '%s%5.3f', ' & ', cp);
        end
        for j=1:numproc
            ml = results.median_length(d,procs(j,1),ih,procs(j,2));
            fprintf(f, '%s%5.3f', ' & ', ml);
        end
        fprintf(f, '%s\n', ' \\');
    end
end

fprintf(f, '%s', '\end{tabular}');

fclose(f);




