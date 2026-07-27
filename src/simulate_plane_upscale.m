function [t_sim_elapsed, t_render_elapsed] = simulate_plane_upscale(varargin)
% SIMULATE_PLANE_UPSCALE  Run a 2D k-Wave plane/line-source simulation and
%                         optionally produce directional plots and an upscaled
%                         rendering video.
%
% USAGE:
%   [t_sim_elapsed, t_render_elapsed] = simulate_plane_upscale(...)
%
% DESCRIPTION:
%   This wrapper builds a 2D k-Wave grid, places a line-like source, runs
%   the time-domain simulation and (optionally) performs a directional gain
%   analysis and/or renders an upscaled video (mp4). The function exposes
%   a number of name-value options to customize domain size, source,
%   frequency, plotting and rendering behaviour.
%
% INPUTS (name-value pairs, defaults shown):
%   'sizex'      : domain size in x (m). Default 0.10
%   'dx'         : grid spacing (m). Default 1e-3
%   'f0'         : center frequency (Hz). Default 40e3
%   'ppw'        : points per wavelength used for dt. Default 10
%   't_end'      : total simulation time (s). Default 600e-6
%   'src_length' : physical length of the planar source (m). Default 10e-3
%   'cx_frac'    : source center x position as fraction of Nx (0..1). Default 0.75
%   'render'     : logical, create mp4 rendering when true. Default false
%   'mp4_filename': output video filename. Default 'kwave_simulation.mp4'
%   'parallel'   : logical, use parfor when creating frames (if available).
%                  Default false
%   'gainAnalysis': logical, perform directional (ring) analysis and plots.
%                   Default false
%   'upscale'    : logical, if true frames are upscaled to target resolution.
%                  Default true
%   'plots'      : cell array or string specifying which plots to produce.
%                  Allowed: 'linear', 'dB', 'interp', 'dbi_abs'. Default {'linear','dB','interp'}
%   'subplotAll' : logical, if true combine requested plots into tiled layout.
%                  Default false
%
% OUTPUTS:
%   t_sim_elapsed    : elapsed time (seconds) for the k-Wave simulation
%   t_render_elapsed : elapsed time (seconds) to render and write the video
%                      (0 if render == false)
%
% NOTES:
%   - Requires the k-Wave MATLAB toolbox (functions like kWaveGrid,
%     kspaceFirstOrder2D, toneBurst).
%   - The directional analysis assumes a 2D line-source geometry and
%     compensates for 1/r geometric spreading when computing relative
%     amplitude in dB.
%   - The video rendering converts frames to an 8-bit RGB movie using a
%     colormap and optional logarithmic mapping for dynamic range.
%
% EXAMPLE:
%   simulate_plane_upscale('sizex',0.2,'dx',2e-3,'f0',50e3,'render',true,'mp4_filename','out.mp4');
%
% AUTHOR:
%   Mikołaj Suchoń — wrapper around k-Wave with convenience plotting/helpers.
%
% ------------------------------------------------------------------------- 

%% --- Input parsing ---
p = inputParser;
addParameter(p, 'sizex', 0.10);
addParameter(p, 'dx', 1e-3);
addParameter(p, 'f0', 40e3);
addParameter(p, 'ppw', 10);
addParameter(p, 't_end', 600e-6);
addParameter(p, 'src_length', 10e-3);
addParameter(p, 'cx_frac', 0.75);
addParameter(p, 'render', false);
addParameter(p, 'mp4_filename', 'kwave_simulation.mp4');
addParameter(p, 'parallel', false);
addParameter(p, 'gainAnalysis', false);
addParameter(p, 'upscale', false);
% default requesting three plot types (linear amplitude, dB, interpolated)
addParameter(p, 'plots', {{'linear','dB','interp'}});
addParameter(p, 'subplotAll', false);
parse(p, varargin{:});
opt = p.Results; 

% When doing gainAnalysis we move the measurement ring closer to source
% (make cx_frac default different to ensure more symmetric sampling).
if opt.gainAnalysis
    opt.cx_frac = 0.5;
end 

% Normalize plot list into a cell array of strings
if ischar(opt.plots)
    opt.plots = {opt.plots};
elseif iscell(opt.plots) && numel(opt.plots) == 1 && iscell(opt.plots{1})
    % allow nested cell form: {'linear','dB'} or {{'linear','dB'}}
    opt.plots = opt.plots{1};
end
if ~iscell(opt.plots)
    opt.plots = {opt.plots};
end 

t_render_elapsed = 0; % default if not rendering 

%% --- Derived parameters ---
sizex = opt.sizex;
dx    = opt.dx;
f0    = opt.f0;
ppw   = opt.ppw;
sizey = sizex;           % square domain
Nx = round(sizex / dx);  % number of grid points in x
Ny = Nx;                 % symmetric grid in y
kgrid = kWaveGrid(Nx, dx, Ny, dx); % build k-Wave grid object 

%% --- Medium setup ---
% Simple homogeneous medium definition (air-like in this example)
medium.sound_speed = 343;     % m/s
medium.density     = 1.225;   % kg/m^3
medium.alpha_coeff = 0;     % absorption coefficient (arbitrary units)
medium.alpha_power = 2;       % frequency power law exponent 

%% --- Time stepping ---
c = medium.sound_speed;
lambda = c / f0;              % wavelength (m)
% choose dt from ppw (points-per-wavelength) and speed c
dt = lambda / (ppw * c);
Nt = round(opt.t_end / dt);   % number of time-steps
kgrid.setTime(Nt, dt);        % set time array in k-Wave grid 

%% --- Source definition ---
% A line-like source (a mask along y at a chosen x index)
source.p_mask = zeros(Nx, Ny);
cx = round(Nx * opt.cx_frac);     % center x index for the source
cy = round(Ny / 2);               % center y index (mid-plane)
half_len = round((opt.src_length / dx) / 2);
yind = (cy - half_len):(cy + half_len);
% ensure indices are inside grid
yind = yind(yind >= 1 & yind <= Ny);
source.p_mask(cx, yind) = 1;      % place the line source mask 

% source temporal waveform — tone burst with number of cycles
num_cycles = 20;
source.p = toneBurst(1/dt, f0, num_cycles); 

%% --- Sensor setup ---
% Record pressure at all grid points
sensor.mask   = true(Nx, Ny);
sensor.record = {'p'}; 

%% --- Simulation (k-Wave) ---
sim_tic = tic;
% pass some options into kspaceFirstOrder2D: single precision, PML, no plots
input_args = {'DataCast','single', 'PMLSize',20, 'PlotPML',false, 'PlotScale',[-0.5,0.5]};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
t_sim_elapsed = toc(sim_tic); 

%% === Directional Gain Analysis (compute p_rms etc) ===
if opt.gainAnalysis
    % sensor_data.p is expected to be [Nx*Ny, Nt]; compute RMS across time
    p_rms = compute_p_rms(sensor_data.p, Nx, Ny);
else
    if ~isempty(opt.plots)
        fprintf('gainAnalysis is false: no directional plots will be produced. Set gainAnalysis=true to enable.\n');
    end
    p_rms = [];
end 

%% === Produce requested directional plots ===
if opt.gainAnalysis && ~isempty(opt.plots)
    % choose a sensible ring radius (in cells) to sample the angular pattern
    radius_cells = choose_ring_radius(lambda, dx, Nx, Ny);
    nAng = 720;
    theta_full = linspace(0,2*pi,nAng+1);
    theta_full(end) = [];  % wrap-around not needed 
% Create interpolant of p_rms on the Nx-by-Ny grid (bilinear interpolation)
F = griddedInterpolant({1:Nx, 1:Ny}, p_rms, 'linear', 'nearest');

% Coordinates (in cell indices) of ring sampling points (xq,yq)
xq = cx + radius_cells * cos(theta_full);
yq = cy + radius_cells * sin(theta_full);

% Keep only points that lie inside the computational grid
inside = xq >= 1 & xq <= Nx & yq >= 1 & yq <= Ny;
if ~all(inside)
    valid_idx = find(inside);
    if isempty(valid_idx)
        error('Sampling ring lies outside grid. Choose smaller radius or increase domain.');
    end
    xq = xq(valid_idx);
    yq = yq(valid_idx);
    theta = theta_full(valid_idx);
else
    theta = theta_full;
end

% Sample interpolated RMS pressure on the ring
meas = F(xq, yq);

% Compute radial distances in meters from center for each sample
r_m = sqrt(((xq - cx) * dx).^2 + ((yq - cy) * dx).^2);

% Normalize measures and compute different representations
linear_norm = abs(meas) / max(abs(meas));    % 0..1 normalized amplitude
gain_lin = abs(meas) .* r_m;                  % compensate for 1/r
gain_lin(gain_lin <= 0) = eps;                % avoid zeros in logs
gain_dB_peak = 20 * log10(gain_lin / max(gain_lin));  % dB normalized to peak = 0
gain_dBi_abs = 20 * log10(gain_lin);                  % absolute dBi-like metric

% Validate requested plot types and normalize to lowercase unique list
requested_plots = unique(cellfun(@lower, opt.plots, 'UniformOutput', false));
allowed = {'linear','db','interp','dbi_abs'};
for i = 1:numel(requested_plots)
    if ~ismember(requested_plots{i}, allowed)
        error('Unknown plot type ''%s''. Allowed: %s', requested_plots{i}, strjoin(allowed,', '));
    end
end

% If user wants all plots in a tiled subplot layout, compute tiling
if opt.subplotAll
    nTiles = 0;
    for pt = requested_plots
        % 'interp' produces two tiles (polar + cartesian)
        if strcmp(pt{1}, 'interp')
            nTiles = nTiles + 2;
        else
            nTiles = nTiles + 1;
        end
    end
    tfig = figure('Name','Directional Patterns (subplots)');
    nrows = max(1, floor(sqrt(nTiles)));
    ncols = ceil(nTiles / nrows);
    tlayout = tiledlayout(nrows, ncols);
    tile_idx = 1;
end

% Loop through requested plot types and render them
for i = 1:numel(requested_plots)
    typ = requested_plots{i};
    switch typ
        case 'linear'
            explain_linear();
            if opt.subplotAll
                nexttile(tlayout, tile_idx);
                polarplot(theta, linear_norm, 'LineWidth', 1.5);
                ax = gca;
                tile_idx = tile_idx + 1;
                title(ax, 'Normalized amplitude (linear) — 0..1, peak at 1');
                set_plot_annotation(ax,...
                    'This polar plot shows the normalized RMS pressure measured on a ring around the source. Values are normalized to the maximum amplitude (peak = 1). Use this to visualize lobes and relative sidelobes in a linear amplitude scale.');
                else
                    figure;
                    polarplot(theta, linear_norm, 'LineWidth', 1.5);
                    title('Normalized amplitude (linear) — 0..1, peak at 1');
                    set_plot_annotation([], ...
'This polar plot shows the normalized RMS pressure measured on a ring around the source. Values are normalized to the maximum amplitude (peak = 1). Use this to visualize lobes and relative sidelobes in a linear amplitude scale.');
                end 
        case 'db'
            explain_dB_peak();
            if opt.subplotAll
                nexttile(tlayout, tile_idx);
                polarplot(theta, gain_dB_peak, 'LineWidth', 1.5);
                ax = gca;
                tile_idx = tile_idx + 1;
                title(ax, 'Relative amplitude (dB), peak = 0 dB');
                try rlim(ax, [-40 0]); catch, end
                set_plot_annotation(ax,...
                    'This polar plot shows amplitude (pressure) on the ring, compensated for 1/r spreading and normalized so the maximum equals 0 dB. This displays lobe widths and sidelobe levels on a dB amplitude scale (20log10).');
                else
                    figure;
                    polarplot(theta, gain_dB_peak, 'LineWidth', 1.5);
                    title('Relative amplitude (dB), peak = 0 dB');
                    try rlim([-40 0]); catch, end
                    set_plot_annotation([], ...
'This polar plot shows amplitude (pressure) on the ring, compensated for 1/r spreading and normalized so the maximum equals 0 dB. This displays lobe widths and sidelobe levels on a dB amplitude scale (20log10).');
                end 
        case 'interp'
            % Interpolated representation + Directivity Index (DI)
            explain_interp_DI();
            if opt.subplotAll
                nexttile(tlayout, tile_idx);
                polarplot(theta, gain_dB_peak, 'LineWidth', 1.25);
                ax1 = gca;
                tile_idx = tile_idx + 1;
                title(ax1, 'Interpolated pattern (polar, peak=0 dB)');
                try rlim(ax1, [-40 0]); catch, end
                set_plot_annotation(ax1,...
                    'Interpolated ring sampling (bilinear) gives a smoother representation of the angular pattern. Compensation for geometric 1/r spreading was applied before normalization.');
                else
                    figure;
                    polarplot(theta, gain_dB_peak, 'LineWidth', 1.25);
                    title('Interpolated pattern (polar, peak = 0 dB)');
                    try rlim([-40 0]); catch, end
                    set_plot_annotation([], ...
'Interpolated ring sampling (bilinear) gives a smoother representation of the angular pattern. Compensation for geometric 1/r spreading was applied before normalization.');
                end 
            % Also provide a Cartesian angle vs dB plot
            angle_deg = rad2deg(theta);
            if opt.subplotAll
                nexttile(tlayout, tile_idx);
                plot(angle_deg, gain_dB_peak, '-k', 'LineWidth', 1.2);
                ax2 = gca;
                tile_idx = tile_idx + 1;
                xlabel(ax2, 'Angle (deg)'); ylabel(ax2, 'Relative amplitude (dB)');
                xlim(ax2, [0 360]); grid(ax2, 'on');
                title(ax2, 'Interpolated pattern (angle vs dB)');
                set_plot_annotation(ax2,...
                    'Cartesian plot makes side-lobe amplitudes and exact angular widths easy to read.');
                else
                    figure;
                    plot(angle_deg, gain_dB_peak, '-k', 'LineWidth', 1.2);
                    xlabel('Angle (deg)'); ylabel('Relative amplitude (dB)');
                    xlim([0 360]); grid on;
                    title('Interpolated pattern (angle vs dB)');
                    set_plot_annotation([], ...
'Cartesian plot makes side-lobe amplitudes and exact angular widths easy to read.');
                end 
            % Compute and print Directivity Index for the 2D case
            DI_dB = compute_DI_2D(gain_lin, theta);
            fprintf('Directivity Index (2D normalization, line-source) = %.2f dB\n', DI_dB);
            fprintf('Physics explanation (DI): DI quantifies how directional the source is. For a line-source (2D) the integral factor is 2π; DI = 10 log10[ (2π * P_max) / ∫ P(θ) dθ ] where P(θ) ∝ (amplitude(θ))^2. Higher DI → more energy concentrated in main lobe.\n');
            
        case 'dbi_abs'
            explain_dBi_abs();
            if opt.subplotAll
                nexttile(tlayout, tile_idx);
                polarplot(theta, gain_dBi_abs, 'LineWidth', 1.5);
                ax = gca;
                tile_idx = tile_idx + 1;
                title(ax, 'Absolute dBi vs 1/r reference (20*log10[p(θ)*r(θ)])');
                set_plot_annotation(ax,...
                    'Absolute dBi (amplitude-based) computed as 20log10( p_rms(θ) * r(θ) ). This compares the measured line-source pressure to a 1/r geometric reference appropriate to a point source at the same location. Note: this is not normalized; offsets depend on absolute amplitudes and measurement radius.');
                else
                    figure;
                    polarplot(theta, gain_dBi_abs, 'LineWidth', 1.5);
                    title('Absolute dBi vs 1/r reference (20log10[p(θ)r(θ)])');
                    set_plot_annotation([],...
                        'Absolute dBi (amplitude-based) computed as 20log10( p_rms(θ) * r(θ) ). This compares the measured line-source pressure to a 1/r geometric reference appropriate to a point source at the same location. Note: this is not normalized; offsets depend on absolute amplitudes and measurement radius.');
                end
        end
    end
end 

%% --- Rendering (video) ---
% If rendering not requested, return now
if ~opt.render
    t_render_elapsed = 0;
    fprintf('Simulation completed (no rendering requested).\n');
    return;
end 

fprintf('Rendering video (upscale = %d)...\n', double(opt.upscale));
[t_render_elapsed] = render_video_frames(sensor_data, Nx, Ny, Nt, opt);
fprintf('Simulation time: %.2f s, Rendering time: %.2f s\n', t_sim_elapsed, t_render_elapsed);
end 

%% ----------------------------- Helper functions ----------------------------
function p_rms = compute_p_rms(pvec, Nx, Ny)
    % compute RMS across time for each spatial location
    % pvec: [Nx*Ny, Nt] matrix (time along 2nd axis)
    p_rms_vec = sqrt(mean(pvec.^2, 2));
    p_rms = reshape(p_rms_vec, [Nx, Ny]);
end 

function r_cells = choose_ring_radius(lambda, dx, Nx, Ny)
    % choose a reasonable sampling ring radius (in grid cells)
    % prefer a radius several wavelengths away but not too close to the domain edge
    min_r_m = max(5 * lambda, 10 * dx);         % at least several wavelengths or several cells
    r_cells = round(min_r_m / dx);
    if r_cells < 3
        % fallback if domain/wavelength make min_r small: choose a fraction of domain
        r_cells = round(0.4 * Nx);
    end
    % also ensure we don't sample too close to domain edges
    r_cells = min([r_cells, floor(0.45 * min(Nx, Ny))]);
end 

function DI_dB = compute_DI_2D(gain_lin, theta)
    % compute Directivity Index (DI) for the 2D line-source case
    % gain_lin is linear amplitude after 1/r compensation (not squared)
    P_theta = gain_lin .^ 2;             % power-like quantity per angle
    P_int = trapz(theta, P_theta);       % integral over angle
    P_max = max(P_theta);                % peak power
    D = (2*pi) * P_max / P_int;          % directivity (2D normalization factor)
    DI_dB = 10 * log10(D);
end 

function out = axifneeded(isSubplot)
    % helper to return gca if requested, otherwise an empty array
    if isSubplot
        out = gca;
    else
        out = [];
    end
end 

function set_plot_annotation(axHandle, textstr)
    % place textual annotation either inside an axis or as a textbox below
    if isempty(axHandle)
        try
            annotation('textbox',[0.13,0.01,0.74,0.08],'String',textstr,'FitBoxToText','on','Interpreter','none','FontSize',9);
        catch
            % annotation may fail in some contexts — ignore silently
        end
    else
        try
            axes(axHandle);
            xlimv = get(axHandle,'XLim'); ylimv = get(axHandle,'YLim');
            tx = text(0.02*(xlimv(2)-xlimv(1))+xlimv(1), 0.02*(ylimv(2)-ylimv(1))+ylimv(1), textstr, 'FontSize', 8, 'Interpreter','none');
            set(tx,'Clipping','off');
        catch
            % ignore any plotting errors
        end
    end
end 

function xlabel_text_description(axHandle)
    % placeholder: could set axis labels consistently if desired
end 

function [t_render_elapsed] = render_video_frames(sensor_data, Nx, Ny, Nt, opt)
    % Render sensor_data.p into a color video and write mp4.
    % sensor_data.p expected [Nx*Ny, Nt]
    t_render_elapsed = 0;
    cmap = parula(256);
    fps = 25; 
% color scaling based on maximum absolute pressure across all frames
clim = max(abs(sensor_data.p(:)));
if clim == 0, clim = 1; end

% target output resolution (height, width) when upscaling
target_res = [1080 1080];

% use a smooth logarithmic mapping to enhance dynamic range in the video
logscale = true;
log_base = 10;

% Pre-allocate cell array for frames
frames = cell(1, Nt);

if opt.parallel
    % parallel frame generation
    parfor it = 1:Nt
        P = reshape(sensor_data.p(:, it), [Nx, Ny]);
        if logscale
            % signed log mapping normalized to [-1,1]
            Pn = sign(P) .* log1p(log_base * abs(P) / clim) / log1p(log_base);
        else
            Pn = P / clim;
        end
        N = 0.5 + 0.5 * Pn;              % map to 0..1
        N = min(max(N, 0), 1);
        idx = max(1, round(N * 255));   % indices into colormap (1..256)
        frame = im2uint8(ind2rgb(idx, cmap));
        if opt.upscale
            frame = imresize(frame, target_res, 'bicubic');
        end
        frames{it} = frame;
    end
else
    % sequential frame generation
    for it = 1:Nt
        P = reshape(sensor_data.p(:, it), [Nx, Ny]);
        if logscale
            Pn = sign(P) .* log1p(log_base * abs(P) / clim) / log1p(log_base);
        else
            Pn = P / clim;
        end
        N = 0.5 + 0.5 * Pn;
        N = min(max(N, 0), 1);
        idx = max(1, round(N * 255));
        frame = im2uint8(ind2rgb(idx, cmap));
        if opt.upscale
            frame = imresize(frame, target_res, 'bicubic');
        end
        frames{it} = frame;
    end
end

% Write frames out to an MP4 using VideoWriter
render_tic = tic;
vw = VideoWriter(opt.mp4_filename, 'MPEG-4');
vw.FrameRate = fps;
open(vw);
for it = 1:Nt
    writeVideo(vw, frames{it});
end
close(vw);
t_render_elapsed = toc(render_tic);

end 

%% ----------------------------- Explanations printed to command window ---------------
function explain_linear()
    fprintf('\n--- Plot: Normalized linear amplitude (polar) ---\n');
    fprintf(['What it shows:\n' ...
'  - RMS pressure on a ring around the source, normalized to the peak (0..1).\n\n']);
end 

function explain_dB_peak()
    fprintf('\n--- Plot: Relative amplitude (dB), peak = 0 dB (polar) ---\n');
    fprintf(['What it shows:\n' ...
'  - Same underlying data as ''linear'', but compensated for 1/r geometric\n' ...
'    spreading and converted to dB; normalized so the maximum is 0 dB.\n\n']);
end 

function explain_interp_DI()
    fprintf('\n--- Plot: Interpolated angular pattern + DI ---\n');
    fprintf(['What it shows:\n' ...
'  - Measurements interpolated smoothly along a ring (bilinear interpolation),\n' ...
'    compensated for 1/r, normalized and shown in dB (polar + Cartesian).\n\n']);
end 

function explain_dBi_abs()
    fprintf('\n--- Plot: Absolute dBi relative to 1/r point-source reference ---\n');
    fprintf(['What it shows:\n' ...
'  - 20log10( p_rms(θ) * r(θ) ). Not normalized — absolute values depend on measurement radius and amplitude.\n\n']);
end
% End of file 
