function t_sim_elapsed = simulate_plane_3d(sizex, dx, f0, mat_filename)
% SIMULATE_PLANE_3D  Run a 3D k-Wave simulation of sound propagation in air
% with a circular piston source in the x–z plane.
%
% Usage:
%   t_sim = simulate_plane_3d()                           % uses defaults
%   t_sim = simulate_plane_3d(sizex, dx, f0, mat_filename)
%
% Defaults:
%   sizex        = 0.10        % [m]
%   dx           = 1e-3        % [m]
%   f0           = 40e3        % [Hz]
%   mat_filename = 'kwave_3d_data.mat'
%
% Output:
%   t_sim_elapsed — total simulation runtime [s]
%
% Requires: k-Wave Toolbox

%— Defaults —%
if nargin < 1 || isempty(sizex), sizex = 0.10; end
if nargin < 2 || isempty(dx), dx = 1e-3; end
if nargin < 3 || isempty(f0), f0 = 40e3; end
if nargin < 4 || isempty(mat_filename), mat_filename = 'kwave_3d_data.mat'; end

%— Grid setup —%
sizey = sizex;
sizez = sizex;
Nx = round(sizex / dx);
Ny = Nx;
Nz = Nx;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

%— Medium (air) —%
medium.sound_speed = 343;   % [m/s]
medium.density     = 1.225; % [kg/m^3]
medium.alpha_coeff = 0.1;   % dB/(MHz^y cm)
medium.alpha_power = 2;

%— Time parameters —%
ppw = 10;                    % points per wavelength
c = medium.sound_speed;
lambda = c / f0;
dt = lambda / (ppw * c);
t_end = 600e-6;              % total time [s]
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

%— Source: circular piston in x–z plane (normal to +y) —%
source.p_mask = zeros(Nx, Ny, Nz);

% Center the source in the middle of x and z, near one side in y
cx = round(Nx / 2);
cy = round(Ny / 4);     % place near one boundary
cz = round(Nz / 2);

radius = round((10e-3) / (2 * dx)); % 1 cm diameter -> 5 mm radius
[xg, yg, zg] = ndgrid(1:Nx, 1:Ny, 1:Nz);

% Create a circular mask in the x–z plane at fixed y = cy
circle_mask = ((xg - cx).^2 + (zg - cz).^2) <= radius^2;
source.p_mask(:, cy, :) = circle_mask(:, cy, :);

% Tone burst excitation
num_cycles = 10;
source.p = toneBurst(1/dt, f0, num_cycles);

%— Sensor: full domain —%
sensor.mask = true(Nx, Ny, Nz);
sensor.record = {'p'};

%— Run simulation —%
disp('Running 3D wave simulation...');
sim_tic = tic;
input_args = {'DataCast', 'single', 'PMLSize', 20, 'PlotPML', false};
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
t_sim_elapsed = toc(sim_tic);

%— Save results —%
disp(['Saving data to ', mat_filename]);
save(mat_filename, 'sensor_data', 'kgrid', 'medium', 'source', 'f0', 'dx', 'sizex', 't_sim_elapsed', '-v7.3');

disp(['Simulation complete. Time elapsed: ', num2str(t_sim_elapsed, '%.2f'), ' s']);
end
