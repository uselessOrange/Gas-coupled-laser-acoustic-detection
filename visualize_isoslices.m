function visualize_isoslices(mat_filename, slice_pos, t_index)
% VISUALIZE_ISOSLICES  Show an isometric slice view of 3D k-Wave data
%
% Usage:
%   visualize_isoslices('kwave_3d_data.mat')
%   visualize_isoslices('kwave_3d_data.mat', [x_slice, y_slice, z_slice])
%   visualize_isoslices('kwave_3d_data.mat', [x_slice, y_slice, z_slice], t_index)
%   visualize_isoslices('kwave_3d_data.mat', [0.05, 0.05, 0.05],100);
%
% Inputs:
%   mat_filename — name of .mat file saved by simulate_plane_3d()
%   slice_pos    — 3-element vector [x y z] in meters (optional)
%   t_index      — time index to visualize (optional, default = end)
%
% The function visualizes the 3D pressure field amplitude using slice planes.

% Load data
data = load(mat_filename);
p = data.sensor_data.p;     % [N_sensor_points x N_time_steps]
kgrid = data.kgrid;
dx = data.dx;

% Extract grid info
Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;
x = (0:Nx-1) * dx;
y = (0:Ny-1) * dx;
z = (0:Nz-1) * dx;

% Pick time index
if nargin < 3 || isempty(t_index)
    t_index = size(p, 2);  % final snapshot
    disp(t_index)
end

% Reshape flattened data back to 3D
p3d = reshape(p(:, t_index), [Nx, Ny, Nz]);

% Compute amplitude (absolute value)
p_amp = abs(p3d);

% Default slice positions at domain center
if nargin < 2 || isempty(slice_pos)
    slice_pos = [x(round(Nx/2)), y(round(Ny/2)), z(round(Nz/2))];
end
x_slice = slice_pos(1);
y_slice = slice_pos(2);
z_slice = slice_pos(3);

% Plot slices
figure('Color', 'w');
xslice = x_slice;
yslice = y_slice;
zslice = z_slice;
h = slice(x, y, z, p_amp, xslice, yslice, zslice);
set(h, 'EdgeColor', 'none', 'FaceColor', 'interp');

colormap(jet);
colorbar;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
title(sprintf('3D Pressure Amplitude (t = %d / %d)', t_index, size(p, 2)));
axis equal tight;
view(3);
camlight;
lighting gouraud;

end
