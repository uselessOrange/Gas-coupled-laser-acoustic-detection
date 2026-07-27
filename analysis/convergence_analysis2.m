function [results, stats] = convergence_analysis2(sizex_vals, dx_vals, f0_vals, render_video)
% CONVERGENCE_ANALYSIS  Sweep dx, sizex, f0 and measure simulation time
% This version runs all simulations in a single parfor loop.

% ----- defaults -----
if nargin < 1 || isempty(sizex_vals), sizex_vals = [0.05, 0.10, 0.20]; end
if nargin < 2 || isempty(dx_vals),    dx_vals    = [2e-3, 1e-3, 5e-4]; end
if nargin < 3 || isempty(f0_vals),    f0_vals    = [5e4, 1e5, 2e5]; end
if nargin < 4 || isempty(render_video), render_video = false; end

% constants (must match simulate_plane)
c = 343;            % sound speed (m/s)
ppw = 20;           % points per wavelength
t_end = 300e-6;

% defaults used for single-variable sweeps
default_sizex = 0.10;
default_dx    = 1e-3;
default_f0    = 1e5;

% confirm simulate_plane exists
if exist('simulate_plane','file') ~= 2
    error('simulate_plane.m not found on MATLAB path. Add it before running this analysis.');
end

% Warm-up (single-threaded) to reduce first-call overhead
fprintf('Warm-up run (no render)...\n');
try
    simulate_plane('sizex', default_sizex, 'dx', default_dx, 'f0', default_f0, 'render', false);
catch ME
    warning('Warm-up simulate failed');
end

% If rendering is requested, fall back to serial loops (rendering in workers is unreliable)
if render_video
    fprintf('render_video == true: running serial execution (rendering cannot reliably run on workers).\n');
    use_parallel = false;
else
    use_parallel = true;
end

% Start parallel pool if requested and not already open
if use_parallel
    if isempty(gcp('nocreate'))
        try
            parpool; % use default profile and number of workers
        catch poolErr
            warning('Could not start parallel pool. Falling back to serial execution.\n');
            use_parallel = false;
        end
    end
end

% Prepare indexing for combined tasks
n_dx = numel(dx_vals);
n_sz = numel(sizex_vals);
n_f  = numel(f0_vals);
total_tasks = n_dx + n_sz + n_f;

% Preallocate combined result arrays
t_all = NaN(total_tasks,1);
err_all = cell(total_tasks,1);
Ncells_all = zeros(total_tasks,1);
Nt_all = zeros(total_tasks,1);
totSamples_all = zeros(total_tasks,1);
task_kind = zeros(total_tasks,1); % 1=dx,2=sizex,3=f0
param_val = zeros(total_tasks,1); % dx or sizex or f0, for convenience

% Fill param arrays (broadcast to workers)
for k = 1:total_tasks
    if k <= n_dx
        task_kind(k) = 1;
        param_val(k) = dx_vals(k);
    elseif k <= n_dx + n_sz
        task_kind(k) = 2;
        param_val(k) = sizex_vals(k - n_dx);
    else
        task_kind(k) = 3;
        param_val(k) = f0_vals(k - n_dx - n_sz);
    end
end

% Choose whether workers should render (never render on workers)
worker_render_flag = false;

% Combined loop (parfor when use_parallel, otherwise for)
if use_parallel
    parfor kk = 1:total_tasks
        try
            % determine which sweep and set parameters
            if task_kind(kk) == 1
                dx = param_val(kk);
                sizex = default_sizex;
                f0 = default_f0;
            elseif task_kind(kk) == 2
                sizex = param_val(kk);
                dx = default_dx;
                f0 = default_f0;
            else
                f0 = param_val(kk);
                sizex = default_sizex;
                dx = default_dx;
            end
            % run simulation on worker (no rendering)
            t_sim_local = simulate_plane('sizex', sizex, 'dx', dx, 'f0', f0, 'render', worker_render_flag);
            t_all(kk) = t_sim_local;
            err_all{kk} = '';

            % derived quantities
            Nx = round(sizex / dx);
            Ny = Nx;
            Ncells_all(kk) = Nx * Ny;
            lambda = c / f0;
            dt = lambda / (ppw * c);       % dt = (lambda/ppw)/c
            Nt = round(t_end / dt);
            Nt_all(kk) = Nt;
            totSamples_all(kk) = Ncells_all(kk) * Nt;

        catch ME
            t_all(kk) = NaN;
            err_all{kk} = ME.message;
            % still compute approximate derived quantities where possible
            try
                if task_kind(kk) == 1
                    dx = param_val(kk);
                    sizex = default_sizex;
                    f0 = default_f0;
                elseif task_kind(kk) == 2
                    sizex = param_val(kk);
                    dx = default_dx;
                    f0 = default_f0;
                else
                    f0 = param_val(kk);
                    sizex = default_sizex;
                    dx = default_dx;
                end
                Nx = round(sizex / dx);
                Ny = Nx;
                Ncells_all(kk) = Nx * Ny;
                lambda = c / f0;
                dt = lambda / (ppw * c);
                Nt = round(t_end / dt);
                Nt_all(kk) = Nt;
                totSamples_all(kk) = Ncells_all(kk) * Nt;
            catch
                % if even that fails, leave zeros/NaNs
            end
        end
    end

else
    % serial execution (prints progress and allows rendering if requested)
    for kk = 1:total_tasks
        if task_kind(kk) == 1
            dx = param_val(kk);
            sizex = default_sizex;
            f0 = default_f0;
            fprintf('Running dx sweep %d/%d: dx=%.4g m\n', kk, total_tasks, dx);
        elseif task_kind(kk) == 2
            sizex = param_val(kk);
            dx = default_dx;
            f0 = default_f0;
            fprintf('Running size sweep %d/%d: sizex=%.3f m\n', kk - n_dx, total_tasks, sizex);
        else
            f0 = param_val(kk);
            sizex = default_sizex;
            dx = default_dx;
            fprintf('Running f0 sweep %d/%d: f0=%.3g Hz\n', kk - n_dx - n_sz, total_tasks, f0);
        end
        try
            t_sim_local = simulate_plane('sizex', sizex, 'dx', dx, 'f0', f0, 'render', render_video);
            t_all(kk) = t_sim_local;
            err_all{kk} = '';
        catch ME
            t_all(kk) = NaN;
            err_all{kk} = ME.message;
            warning('Task %d failed: %s', kk, ME.message);
        end

        % derived quantities
        Nx = round(sizex / dx);
        Ny = Nx;
        Ncells_all(kk) = Nx * Ny;
        lambda = c / f0;
        dt = lambda / (ppw * c);       % dt = (lambda/ppw)/c
        Nt = round(t_end / dt);
        Nt_all(kk) = Nt;
        totSamples_all(kk) = Ncells_all(kk) * Nt;
    end

end

% Now split combined results back into per-sweep arrays
dx_t = t_all(1:n_dx);
dx_err = err_all(1:n_dx);
dx_Ncells = Ncells_all(1:n_dx);
dx_Nt = Nt_all(1:n_dx);
dx_totSamples = totSamples_all(1:n_dx);

sz_t = t_all(n_dx+1 : n_dx+n_sz);
sz_err = err_all(n_dx+1 : n_dx+n_sz);
sz_Ncells = Ncells_all(n_dx+1 : n_dx+n_sz);
sz_Nt = Nt_all(n_dx+1 : n_dx+n_sz);
sz_totSamples = totSamples_all(n_dx+1 : n_dx+n_sz);

f_t = t_all(n_dx+n_sz+1 : end);
f_err = err_all(n_dx+n_sz+1 : end);
f_Ncells = Ncells_all(n_dx+n_sz+1 : end);
f_Nt = Nt_all(n_dx+n_sz+1 : end);
f_totSamples = totSamples_all(n_dx+n_sz+1 : end);

% Build results tables (same layout as before)
results.dx_sweep = table(dx_vals(:), dx_t, dx_Ncells, dx_Nt, dx_totSamples, dx_err, ...
    'VariableNames', {'dx','t_sim','Ncells','Nt','TotalSamples','Error'});

results.size_sweep = table(sizex_vals(:), sz_t, sz_Ncells, sz_Nt, sz_totSamples, sz_err, ...
    'VariableNames', {'sizex','t_sim','Ncells','Nt','TotalSamples','Error'});

results.f0_sweep = table(f0_vals(:), f_t, f_Ncells, f_Nt, f_totSamples, f_err, ...
    'VariableNames', {'f0','t_sim','Ncells','Nt','TotalSamples','Error'});

% -------------------------
% Combine all runs for regression (exclude NaNs)
% -------------------------
all_t = [dx_t; sz_t; f_t];
all_Ncells = [dx_Ncells; sz_Ncells; f_Ncells];
all_Nt = [dx_Nt; sz_Nt; f_Nt];
all_total = [dx_totSamples; sz_totSamples; f_totSamples];
valid = ~isnan(all_t) & all_total>0;

if sum(valid) < 2
    warning('Not enough valid samples to fit regression models.');
    stats.lm = [];
    stats.lm_total = [];
    stats.per_sample = struct('mean',NaN,'median',NaN,'std',NaN);
    stats.per_sample_time = NaN(size(all_t));
else
    % Fit linear model: t_sim = aNcells + bNt + c
    tbl = table(all_t(valid), all_Ncells(valid), all_Nt(valid), 'VariableNames', {'t_sim','Ncells','Nt'});
    lm = fitlm(tbl, 't_sim ~ Ncells + Nt');
    % Also compute simple regression vs total samples
    lm_total = fitlm(all_total(valid), all_t(valid));

    % per-sample time
    per_sample_time = all_t(valid) ./ all_total(valid); % seconds per sample
    per_sample_stats.mean = mean(per_sample_time);
    per_sample_stats.median = median(per_sample_time);
    per_sample_stats.std = std(per_sample_time);

    stats.lm = lm;
    stats.lm_total = lm_total;
    stats.per_sample = per_sample_stats;
    stats.per_sample_time = per_sample_time;

end

% Keep summary table (including NaNs and per-iteration errors)
results.summary_all = table(all_t, all_Ncells, all_Nt, all_total, 'VariableNames', ...
    {'t_sim','Ncells','Nt','TotalSamples'});

% also attach per-sweep error messages
stats.errors.dx = dx_err;
stats.errors.sizex = sz_err;
stats.errors.f0 = f_err;

% Plotting (works with NaNs)
figure('Name','Convergence vs time','NumberTitle','off');
% dx vs time
subplot(2,2,1);
loglog(dx_vals, dx_t, 'o-','LineWidth',1.4);
xlabel('dx (m)');
ylabel('t_{sim} (s)');
title('dx vs t_{sim}');
grid on;

% sizex vs time
subplot(2,2,2);
loglog(sizex_vals, sz_t, 's-','LineWidth',1.4);
xlabel('sizex (m)');
ylabel('t_{sim} (s)');
title('sizex vs t_{sim}');
grid on;

% f0 vs time
subplot(2,2,3);
loglog(f0_vals, f_t, 'd-','LineWidth',1.4);
xlabel('f_0 (Hz)');
ylabel('t_{sim} (s)');
title('f_0 vs t_{sim}');
grid on;

% Annotated summary (regression)
subplot(2,2,4);
axis off;
if isfield(stats,'lm') && ~isempty(stats.lm)
    % extract coefficients robustly by row name (works for fitlm)
    coefTable = stats.lm.Coefficients;
    rowNames = coefTable.Properties.RowNames;
    % a = coefficient for Ncells
    idx_a = find(strcmp(rowNames,'Ncells'));
    idx_b = find(strcmp(rowNames,'Nt'));
    idx_c = find(strcmp(rowNames,'(Intercept)'));
    if isempty(idx_a), a = NaN; else a = coefTable.Estimate(idx_a); end
    if isempty(idx_b), b = NaN; else b = coefTable.Estimate(idx_b); end
    if isempty(idx_c), c0 = NaN; else c0 = coefTable.Estimate(idx_c); end
    R2 = stats.lm.Rsquared.Ordinary;
    slope_total = stats.lm_total.Coefficients.Estimate(2);
    per_mean = stats.per_sample.mean;
    per_med = stats.per_sample.median;
    per_std = stats.per_sample.std;
    txt = {
        sprintf('Linear model: t_sim = aNcells + bNt + c')
        sprintf('a = %.3e  (s per cell)', a)
        sprintf('b = %.3e  (s per time-step)', b)
        sprintf('c = %.3e  (s)         ', c0)
        sprintf('R^2 = %.3f', R2)
        ''
        sprintf('Total-sample regression slope: %.3e s/sample', slope_total)
        ''
        sprintf('Per-sample time (mean/median/std): %.3e / %.3e / %.3e s', per_mean, per_med, per_std)
        };
else
    txt = {'Not enough valid samples to fit regression models.'};
end
text(0, 1, txt, 'VerticalAlignment','top', 'FontName','FixedWidth');
title('Regression summary');

fprintf('Analysis complete. See figure and results struct.\n');

end
