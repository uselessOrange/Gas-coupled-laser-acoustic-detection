function visualize_plane_3d(mat_filename)
% VISUALIZE_PLANE_3D  Visualize saved 3D k-Wave simulation data
%
% Usage:
%   visualize_plane_3d('kwave_3d_data.mat')

if nargin < 1 || isempty(mat_filename)
    mat_filename = 'kwave_3d_data.mat';
end

%— Load data —%
load(mat_filename, 'sensor_data', 'kgrid');

%— Extract parameters —%
p = sensor_data.p;
[Nx, Ny, Nz, Nt] = size(p);

%— Visualization setup —%
figure('Color', 'w');
colormap(jet);
clim = max(abs(p(:)));
if clim == 0, clim = 1; end

disp('Press Ctrl+C to stop visualization.');
for it = 1:Nt
    slice_x = round(Nx/2);
    slice_y = round(Ny/2);
    slice_z = round(Nz/2);

    % Extract central slices
    Px = squeeze(p(slice_x, :, :, it));
    Py = squeeze(p(:, slice_y, :, it));
    Pz = squeeze(p(:, :, slice_z, it));

    % Plot slices
    subplot(1,3,1);
    imagesc(Px'); axis image off; title('X Slice'); caxis([-clim clim]);
    subplot(1,3,2);
    imagesc(Py'); axis image off; title('Y Slice'); caxis([-clim clim]);
    subplot(1,3,3);
    imagesc(Pz'); axis image off; title('Z Slice'); caxis([-clim clim]);
    sgtitle(sprintf('Time step %d / %d', it, Nt));
    drawnow;
end
end
