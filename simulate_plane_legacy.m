function [t_sim_elapsed, t_render_elapsed] = simulate_plane_legacy(...
        sizex, dx, f0, render_video, mp4_filename,parallel)
% SIMULATE_PLANE  Run a 2D k-Wave simulation on a square plane
%
% Usage:
%   [t_sim, t_render] = simulate_plane()                        % uses all defaults
%   [t_sim, t_render] = simulate_plane(sizex)                  % override sizex
%   [t_sim, t_render] = simulate_plane(sizex, dx, f0)          % override first 3
%   [t_sim, t_render] = simulate_plane(sizex, dx, f0, render)  % override all
%
% Defaults (from your script):
%   sizex        = 0.10      % [m]
%   dx           = 1e-3      % [m]
%   f0           = 40e3       % [Hz]
%   render_video = true      % produce MP4

%— set defaults if not provided —%
if nargin<1 || isempty(sizex),        sizex = 0.10;   end
if nargin<2 || isempty(dx),           dx    = 1e-3;   end
if nargin<3 || isempty(f0),           f0    = 40e3;    end
if nargin<4 || isempty(render_video), render_video = false; end
if nargin<5 || isempty(mp4_filename), mp4_filename = 'kwave_simulation.mp4'; end
if nargin<6 || isempty(parallel), parallel = false; end
    
t_render_elapsed=0;    
%— derived grid parameters —%
sizey = sizex;
Nx    = round(sizex/dx);
Ny    = Nx;
kgrid = kWaveGrid(Nx, dx, Ny, dx);

%— medium —%
medium.sound_speed = 343;   % m/s
medium.density     = 1.225; % kg/m^3
medium.alpha_coeff = 0.1;   % dB/(MHz^y cm)
medium.alpha_power = 2;

%— time stepping —%
ppw    = 10;                    % points per wavelength
c      = medium.sound_speed;
lambda = c / f0;
dt     = lambda / (ppw * c);
t_end  = 600e-6;                % total time [s]
Nt     = round(t_end / dt);
kgrid.setTime(Nt, dt);

%— source: horizontal 1 cm tone burst —%
source.p_mask = zeros(Nx,Ny);
cx = round(Nx-Nx/4);
cy = round(Ny/2);
half_len = round((10e-3)/dx/2);  % ±5 mm → 1 cm total
yind     = (cy-half_len):(cy+half_len);
source.p_mask(cx, yind) = 1;
num_cycles = 20;
source.p    = toneBurst(1/dt, f0, num_cycles);

%— sensor: record entire field —%
sensor.mask   = true(Nx,Ny);
sensor.record = {'p'};

%— run simulation —%
sim_tic = tic;
input_args = {'DataCast','single', 'PMLSize',20, 'PlotPML',false, 'PlotScale',[-0.5,0.5]};
sensor_data = kspaceFirstOrder2D( ...
    kgrid, medium, source, sensor, input_args{:} );
t_sim_elapsed = toc(sim_tic);

%— if rendering skipped —%
if ~render_video
    t_render_elapsed = 0;
    return
end

%— prepare video rendering —%
fps          = 25;
clim         = max(abs(sensor_data.p(:)));
if clim==0, clim=1; end
cmap = parula(256);

%— precompute RGB frames in parallel —%
fprintf('Generating %d frames in parallel...\n', Nt);
frames = cell(1, Nt);
if parallel

    parfor it = 1:Nt
        P = reshape(sensor_data.p(:,it), [Nx,Ny]);
        N = (P + clim) / (2*clim);
        N = min(max(N,0),1);
        idx = max(1, round(N*255));
        frames{it} = im2uint8(ind2rgb(idx, cmap));
    end
else
    for it = 1:Nt
        P = reshape(sensor_data.p(:,it), [Nx,Ny]);
        N = (P + clim) / (2*clim);
        N = min(max(N,0),1);
        idx = max(1, round(N*255));
        frames{it} = im2uint8(ind2rgb(idx, cmap));
    end
end
%— write MP4 —%
render_tic = tic;
vw = VideoWriter(mp4_filename, 'MPEG-4');
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