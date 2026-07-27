function [t_sim_elapsed, t_render_elapsed] = simulate_plane_upscale_legacy2(varargin)
% SIMULATE_PLANE_UPSCALE  Run a 2D k-Wave simulation with optional upscaling.
%
% Usage:
%   simulate_plane_upscale()                              % defaults
%   simulate_plane_upscale('gainAnalysis',true,'upscale',false)
%
% New / changed inputs:
%   upscale      = true/false   % if true keep previous upscaling behavior
%   gainAnalysis = true/false   % if true enforce cx_frac = 0.5 (source centered)
%
% Returns:
%   t_sim_elapsed, t_render_elapsed 

%% --- Input Parsing ---
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
addParameter(p, 'upscale', true);  % new: whether to upscale frames for rendering
parse(p, varargin{:});
opt = p.Results; 

% If gain analysis requested, force source to be centered in x
if opt.gainAnalysis
    opt.cx_frac = 0.5;
end 

t_render_elapsed = 0; 

%% --- Derived parameters ---
sizex = opt.sizex;
dx    = opt.dx;
f0    = opt.f0;
ppw   = opt.ppw;
sizey = sizex;
Nx = round(sizex / dx);
Ny = Nx;
kgrid = kWaveGrid(Nx, dx, Ny, dx); 

%% --- Medium setup ---
medium.sound_speed = 343;   % m/s
medium.density     = 1.225; % kg/m^3
medium.alpha_coeff = 0.1;   % dB/(MHz^y cm)
medium.alpha_power = 2; 

%% --- Time stepping ---
c = medium.sound_speed;
lambda = c / f0;
dt = lambda / (ppw * c);
Nt = round(opt.t_end / dt);
kgrid.setTime(Nt, dt); 

%% --- Source definition ---
source.p_mask = zeros(Nx, Ny);
cx = round(Nx * opt.cx_frac);
cy = round(Ny / 2);
half_len = round((opt.src_length / dx) / 2);
yind = (cy - half_len):(cy + half_len);
yind = yind(yind >= 1 & yind <= Ny);  % clamp indices
source.p_mask(cx, yind) = 1;
num_cycles = 20;
source.p = toneBurst(1/dt, f0, num_cycles); 

%% --- Sensor setup ---
sensor.mask   = true(Nx, Ny);
sensor.record = {'p'}; 

%% --- Simulation ---
sim_tic = tic;
input_args = {'DataCast','single', 'PMLSize',20, 'PlotPML',false, 'PlotScale',[-0.5,0.5]};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
t_sim_elapsed = toc(sim_tic); 

%% === Directional Gain Analysis: normalized amplitude in dB + DI ===
if opt.gainAnalysis
    % time-RMS pressure at each sensor (reshape to grid)
    p_rms = sqrt(mean(sensor_data.p.^2, 2));
    p_rms = reshape(p_rms, [Nx, Ny]); 
% source center (cx_frac forced to 0.5 if opt.gainAnalysis set earlier)
cx = round(Nx * opt.cx_frac);
cy = round(Ny / 2);

% choose a radius in grid cells — pick far field: >= 5..10 wavelengths
lambda = c / f0;
min_r_m = max(5 * lambda, 10 * dx);           % at least several lambda
radius = round(min_r_m / dx);
% if that is too small for grid, use a fraction of domain
if radius < 3
    radius = round(0.4 * Nx);
end

% angular sampling
nAng = 720;
theta = linspace(0, 2*pi, nAng+1); theta(end) = [];
% continuous (floating) sampling positions in grid indices
xq = cx + (radius * cos(theta));   % row index (x direction)
yq = cy + (radius * sin(theta));   % col index (y direction)

% ensure points fall inside grid; if some outside, shrink radius slightly
inside = xq >= 1 & xq <= Nx & yq >= 1 & yq <= Ny;
if ~all(inside)
    valid_idx = find(inside);
    if isempty(valid_idx)
        error('Sampling ring lies fully outside grid. Choose a smaller radius.');
    end
    xq = xq(valid_idx); yq = yq(valid_idx); theta = theta(valid_idx);
end

% bilinear interpolation using griddedInterpolant for smooth sampling
F = griddedInterpolant({1:Nx, 1:Ny}, p_rms, 'linear', 'nearest');
meas = F(xq, yq);               % measured RMS pressures at ring (amplitude)

% distances in meters from the source center to ring points
r_m = sqrt(((xq - cx) * dx).^2 + ((yq - cy) * dx).^2);

% compensate for 1/r of point-source reference → amplitude comparison
gain_lin = abs(meas) .* r_m;

% avoid zeros, normalize so peak = 0 dB, convert to dB (amplitude → 20*log10)
gain_lin(gain_lin <= 0) = eps;
gain_dB = 20 * log10(gain_lin / max(gain_lin));

% plot: polar with radial limit for readability
figure;
polarplot(theta, gain_dB, 'LineWidth', 1.5);
title('Directional pattern (amplitude, peak = 0 dB) relative to 1/r point source');
try
    rlim([-40 0]);  % show down to -40 dB (adjust as needed)
catch
    % older MATLAB versions may not support rlim for polar plots
end
% also show a Cartesian angular plot for precise sidelobe levels
figure;
plot(rad2deg(theta), gain_dB, '-k', 'LineWidth', 1.2);
xlabel('Angle (deg)'); ylabel('Relative amplitude (dB)');
xlim([0 360]); grid on; title('Directional pattern vs angle');

% Optional: compute a 2D directivity index (DI) for the planar / 2D case.
% NOTE: your simulation is 2D (line source), so the standard 3D factor 4*pi
% becomes 2*pi for 2D angular integration. We compute:
%   D = (2*pi) * P_max / ∫_0^{2π} P(θ) dθ, where P ∝ (gain_lin)^2 is power
P_theta = gain_lin .^ 2;
P_int = trapz(theta, P_theta);        % integral over 0..2π
P_max = max(P_theta);
D = (2*pi) * P_max / P_int;
DI_dB = 10 * log10(D);
fprintf('Directivity Index (2D, DI) = %.2f dB (note: 2D normalization used)\n', DI_dB);

end 
%% --- Rendering ---
if ~opt.render
    t_render_elapsed = 0;
    fprintf('Simulation completed (no rendering).\n');
    return;
end 

fprintf('Rendering video%s...\n', ternary(opt.upscale, ' (Full HD upscale)', '')); 

fps = 25;
clim = max(abs(sensor_data.p(:)));
if clim == 0, clim = 1; end
cmap = parula(256); 

% upscale target resolution (only used if opt.upscale true)
target_res = [1080 1080];   % Full HD (square chosen to match previous code) 

% mild dynamic range compression
logscale = true;
log_base = 10; 

frames = cell(1, Nt); 

% nested helper to convert a pressure snapshot to an RGB frame
    function frame = make_frame(P)
        % P: Nx-by-Ny pressure snapshot
        if logscale
            Pn = sign(P) .* log1p(log_base * abs(P) / clim) / log1p(log_base);
        else
            Pn = P / clim;
        end
        % normalize -1..1 -> 0..1
        N = 0.5 + 0.5 * Pn;
        N = min(max(N, 0), 1);
        idx = max(1, round(N * 255));
        frame = im2uint8(ind2rgb(idx, cmap));
        if opt.upscale
            frame = imresize(frame, target_res, 'bicubic');
        end
    end 

if opt.parallel
    % parfor requires frames to be sliced and nested functions may not be supported
    % in all MATLAB versions inside parfor. We keep it as close as possible, but if
    % nested make_frame fails in your MATLAB version inside parfor, switch to serial.
    parfor it = 1:Nt
        P = reshape(sensor_data.p(:, it), [Nx, Ny]);
        % inline the conversion if nested functions fail in parfor
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
else
    for it = 1:Nt
        P = reshape(sensor_data.p(:, it), [Nx, Ny]);
        frames{it} = make_frame(P);
    end
end 

render_tic = tic;
vw = VideoWriter(opt.mp4_filename, 'MPEG-4');
vw.FrameRate = fps;
open(vw);
for it = 1:Nt
    writeVideo(vw, frames{it});
end
close(vw);
t_render_elapsed = toc(render_tic); 

fprintf('Simulation time: %.2f s, Rendering time: %.2f s\n', ...
    t_sim_elapsed, t_render_elapsed);
end 

% small utility ternary (inline to avoid external dependencies)
function out = ternary(cond, a, b)
if cond, out = a; else out = b; end
end 
