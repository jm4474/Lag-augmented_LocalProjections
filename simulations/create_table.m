clear;

% Create table of coverage and length


%% Settings

% DGP type
dgp_type = 'iid'; % Either 'iid', 'arch', or 'homosk';

% File names
load_filename = strcat('../results/sim_', dgp_type, '.mat'); % Load results from this file
save_filename = strcat('./', dgp_type); % First part of file name for saved figures

% Select DGPs
T_sel = 240;
rhos_sel = [0.0 0.5 0.95 1];

%{
% Select CI procedures
if strcmp(se_type, 'boot')
    proc_names = {'$\text{AR}_d$', '$\text{AR}_b^*$', '$\text{LP}_b$', '$\text{LP}_b^*$'}; % Names of procedures
    procs = [1 1;
             4 2;
             6 4;
             7 4]; % Rows = selected procedures; 1st column = index of procedure; 2nd column = index of CI type
else
    proc_names = {'$\text{AR}_d$', '$\text{AR}_d^*$', '$\text{LP}_d$', '$\text{LP}_d^*$'};
    procs = [1 1;
             4 1;
             6 1;
             7 1];
end
%}
proc_names = {'$\text{LP-LA}_b$', '$\text{LP-LA}$', '$\text{LP}_b$', '$\text{LP}$', '$\text{AR-LA}_b$', '$\text{AR}$'};
procs = [7 4;
         7 1;
         6 4;
         6 1
         4 2
         1 1];
%{
proc_names = {'$\text{AR}_e$', '$\text{LP}_h$', '$\text{LP}_w$', '$\text{LP}_h$'};
procs = [5 2;
         6 1;
         7 1;
         7 4];
%}

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

f = fopen(strcat(save_filename, '.tex'), 'w'); % Open file for writing

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
            color_str = '';
            fprintf(f, '%s%s%5.3f%s', ' & ', '{', cp, '}');
        end
        for j=1:numproc
            ml = results.median_length(d,procs(j,1),ih,procs(j,2));
            color_str = '';
            fprintf(f, '%s%s%s%5.3f%s', ' & ', color_str, '{', ml, '}');
        end
        fprintf(f, '%s\n', ' \\');
    end
end

fprintf(f, '%s', '\end{tabular}');

fclose(f);




