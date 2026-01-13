function [results, stats] = convergence_analysis(sizex_vals, dx_vals, f0_vals, render_video)
% CONVERGENCE_ANALYSIS  Sweep dx, sizex, f0 and measure simulation time (parfor-enabled)
%
% See your original function for usage and defaults. This version uses parfor
% to run the independent sweep iterations in parallel when render_video==false.

% ----- defaults -----
if nargin < 1 || isempty(sizex_vals), sizex_vals = [0.05, 0.10, 0.20]; end
if nargin < 2 || isempty(dx_vals),    dx_vals    = [2e-3, 1e-3, 5e-4]; end
if nargin < 3 || isempty(f0_vals),    f0_vals    = [5e4, 1e5, 2e5]; end
if nargin < 4 || isempty(render_video), render_video = false; end

% constants (must match simulate_plane)
c = 343;            % sound speed
ppwc = 20;
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
    simulate_plane(default_sizex, default_dx, default_f0, false);
catch ME
    warning('Warm-up simulate failed:');
end

% If rendering is requested, fall back to serial loops (rendering in workers is unreliable)
if render_video
    fprintf('render_video == true: running serial loops (rendering cannot reliably run on workers).\n');
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
            warning('Could not start parallel pool \nFalling back to serial execution.');
            use_parallel = false;
        end
    end
end

% -------------------------
% 1) Sweep dx (vary spatial step)
% -------------------------
n_dx = numel(dx_vals);
dx_t = NaN(n_dx,1);            % simulation time
dx_err = cell(n_dx,1);        % error messages
dx_Ncells = zeros(n_dx,1);
dx_Nt = zeros(n_dx,1);
dx_totSamples = zeros(n_dx,1);

if use_parallel
    parfor k=1:n_dx
        dx = dx_vals(k);
        sizex = default_sizex;
        f0 = default_f0;
        % run simulation on worker
        try
            % ensure we do not render inside workers
            t_sim = simulate_plane(sizex, dx, f0, false);
            dx_t(k) = t_sim;
            dx_err{k} = '';
        catch ME
            dx_t(k) = NaN;
            dx_err{k} = ME.message;
        end
        Nx = round(sizex/dx); Ny = Nx;
        dx_Ncells(k) = Nx * Ny;
        lambda = c/f0;
        dt = lambda/(ppwc);
        Nt = round(t_end/dt);
        dx_Nt(k) = Nt;
        dx_totSamples(k) = dx_Ncells(k) * Nt;
    end
else
    for k=1:n_dx
        dx = dx_vals(k);
        sizex = default_sizex;
        f0 = default_f0;
        fprintf('Running dx sweep %d/%d: dx=%.4g m\n', k, n_dx, dx);
        try
            t_sim = simulate_plane(sizex, dx, f0, render_video);
            dx_t(k) = t_sim;
            dx_err{k} = '';
        catch ME
            dx_t(k) = NaN;
            dx_err{k} = ME.message;
            warning('dx sweep iteration %d failed: %s', k, ME.message);
        end
        Nx = round(sizex/dx); Ny = Nx;
        dx_Ncells(k) = Nx * Ny;
        lambda = c/f0;
        dt = lambda/(ppwc);
        Nt = round(t_end/dt);
        dx_Nt(k) = Nt;
        dx_totSamples(k) = dx_Ncells(k) * Nt;
    end
end
results.dx_sweep = table(dx_vals(:), dx_t, dx_Ncells, dx_Nt, dx_totSamples, dx_err, ...
    'VariableNames', {'dx','t_sim','Ncells','Nt','TotalSamples','Error'});

% -------------------------
% 2) Sweep sizex (vary plane size)
% -------------------------
n_sz = numel(sizex_vals);
sz_t = NaN(n_sz,1);
sz_err = cell(n_sz,1);
sz_Ncells = zeros(n_sz,1);
sz_Nt = zeros(n_sz,1);
sz_totSamples = zeros(n_sz,1);

if use_parallel
    parfor k=1:n_sz
        sizex = sizex_vals(k);
        dx = default_dx;
        f0 = default_f0;
        try
            t_sim = simulate_plane(sizex, dx, f0, false);
            sz_t(k) = t_sim;
            sz_err{k} = '';
        catch ME
            sz_t(k) = NaN;
            sz_err{k} = ME.message;
        end
        Nx = round(sizex/dx); Ny = Nx;
        sz_Ncells(k) = Nx * Ny;
        lambda = c/f0;
        dt = lambda/(ppwc);
        Nt = round(t_end/dt);
        sz_Nt(k) = Nt;
        sz_totSamples(k) = sz_Ncells(k) * Nt;
    end
else
    for k=1:n_sz
        sizex = sizex_vals(k);
        dx = default_dx;
        f0 = default_f0;
        fprintf('Running size sweep %d/%d: sizex=%.3f m\n', k, n_sz, sizex);
        try
            t_sim = simulate_plane(sizex, dx, f0, render_video);
            sz_t(k) = t_sim;
            sz_err{k} = '';
        catch ME
            sz_t(k) = NaN;
            sz_err{k} = ME.message;
            warning('size sweep iteration %d failed: %s', k, ME.message);
        end
        Nx = round(sizex/dx); Ny = Nx;
        sz_Ncells(k) = Nx * Ny;
        lambda = c/f0;
        dt = lambda/(ppwc);
        Nt = round(t_end/dt);
        sz_Nt(k) = Nt;
        sz_totSamples(k) = sz_Ncells(k) * Nt;
    end
end
results.size_sweep = table(sizex_vals(:), sz_t, sz_Ncells, sz_Nt, sz_totSamples, sz_err, ...
    'VariableNames', {'sizex','t_sim','Ncells','Nt','TotalSamples','Error'});

% -------------------------
% 3) Sweep f0 (vary frequency)
% -------------------------
n_f = numel(f0_vals);
f_t = NaN(n_f,1);
f_err = cell(n_f,1);
f_Ncells = zeros(n_f,1);
f_Nt = zeros(n_f,1);
f_totSamples = zeros(n_f,1);

if use_parallel
    parfor k=1:n_f
        f0 = f0_vals(k);
        sizex = default_sizex;
        dx = default_dx;
        try
            t_sim = simulate_plane(sizex, dx, f0, false);
            f_t(k) = t_sim;
            f_err{k} = '';
        catch ME
            f_t(k) = NaN;
            f_err{k} = ME.message;
        end
        Nx = round(sizex/dx); Ny = Nx;
        f_Ncells(k) = Nx * Ny;
        lambda = c/f0;
        dt = lambda/(ppwc);
        Nt = round(t_end/dt);
        f_Nt(k) = Nt;
        f_totSamples(k) = f_Ncells(k) * Nt;
    end
else
    for k=1:n_f
        f0 = f0_vals(k);
        sizex = default_sizex;
        dx = default_dx;
        fprintf('Running f0 sweep %d/%d: f0=%.3g Hz\n', k, n_f, f0);
        try
            t_sim = simulate_plane(sizex, dx, f0, render_video);
            f_t(k) = t_sim;
            f_err{k} = '';
        catch ME
            f_t(k) = NaN;
            f_err{k} = ME.message;
            warning('f0 sweep iteration %d failed: %s', k, ME.message);
        end
        Nx = round(sizex/dx); Ny = Nx;
        f_Ncells(k) = Nx * Ny;
        lambda = c/f0;
        dt = lambda/(ppwc);
        Nt = round(t_end/dt);
        f_Nt(k) = Nt;
        f_totSamples(k) = f_Ncells(k) * Nt;
    end
end
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
    a = stats.lm.Coefficients.Estimate('Ncells');
    b = stats.lm.Coefficients.Estimate('Nt');
    c0 = stats.lm.Coefficients.Estimate('(Intercept)');
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