function [t_sim_elapsed, t_render_elapsed] = simulate_plane(varargin)
% SIMULATE_PLANE  Run a 2D k-Wave simulation on a square plane with flexible options.
%
% Usage:
%   simulate_plane()                              % all defaults
%   simulate_plane('sizex', 0.15, 'f0', 25e3)     % override some parameters
%   simulate_plane('render', false, 'parallel', true)
%
% Parameters (defaults):
%   sizex        = 0.10       % [m]
%   dx           = 1e-3       % [m]
%   f0           = 40e3       % [Hz]
%   ppw          = 10         % points per wavelength
%   src_length   = 10e-3      % total source length [m]
%   cx_frac      = 0.75       % source x-position as fraction of Nx (1 - 1/4)
%   render       = true       % whether to render video
%   mp4_filename = 'kwave_simulation.mp4'
%   parallel     = false      % parallel rendering
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
addParameter(p, 'render', true);
addParameter(p, 'mp4_filename', 'kwave_simulation.mp4');
addParameter(p, 'parallel', false);
parse(p, varargin{:});
opt = p.Results;

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
%t_end = 600e-6; % total time [s]
Nt = round(opt.t_end / dt);
kgrid.setTime(Nt, dt);

%% --- Source definition ---
source.p_mask = zeros(Nx, Ny);
cx = round(Nx * opt.cx_frac);
cy = round(Ny / 2);
half_len = round((opt.src_length / dx) / 2);
yind = (cy - half_len):(cy + half_len);
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

%% --- Rendering ---
if ~opt.render
    t_render_elapsed = 0;
    fprintf('Simulation completed (no rendering).\n');
    return;
end

fprintf('Rendering video...\n');
fps = 25;
clim = max(abs(sensor_data.p(:)));
if clim == 0, clim = 1; end
cmap = parula(256);

frames = cell(1, Nt);
if opt.parallel
    parfor it = 1:Nt
        P = reshape(sensor_data.p(:, it), [Nx, Ny]);
        N = (P + clim) / (2 * clim);
        N = min(max(N, 0), 1);
        idx = max(1, round(N * 255));
        frames{it} = im2uint8(ind2rgb(idx, cmap));
    end
else
    for it = 1:Nt
        P = reshape(sensor_data.p(:, it), [Nx, Ny]);
        N = (P + clim) / (2 * clim);
        N = min(max(N, 0), 1);
        idx = max(1, round(N * 255));
        frames{it} = im2uint8(ind2rgb(idx, cmap));
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
