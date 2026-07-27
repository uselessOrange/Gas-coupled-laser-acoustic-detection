clear
tic


% -------------------- GRID / MEDIUM / SOURCE (same as before) ------------
% Nx = 256; Ny = 256; %assume (256x256)/5 mm plane
% dx = (1e-3)/5; dy = (1e-3)/5; % step is 0.05mm
% kgrid = kWaveGrid(Nx, dx, Ny, dy);

sizex= 0.1; sizey=sizex; %10cm


dx = 1e-3; dy = dx; % step is 0.1mm
Nx = sizex/dx; Ny = Nx; 

kgrid = kWaveGrid(Nx, dx, Ny, dy);


medium.sound_speed = 343;
medium.density = 1.225;
medium.alpha_coeff = 0.1;
medium.alpha_power = 2;

f0 = 1e5;                % 100 kHz
ppw = 20;
c = medium.sound_speed;
lambda = c / f0;
dt = lambda / (ppw * c);
t_end = 300e-6;
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% -------------------- SOURCE: horizontal line, 1 cm long -----------------
source.p_mask = zeros(Nx, Ny);

center_x = round(Nx/2);
center_y = round(Ny/2);
half_length = round((10e-3 / dy) / 2); % half of 1 cm = 5 mm = 5 grid points

% Create a horizontal line centered in the grid
y_indices = (center_y - half_length):(center_y + half_length);
source.p_mask(center_x, y_indices) = 1;

num_cycles = 10;
source.p = toneBurst(1/kgrid.dt, f0, num_cycles);

% -------------------- SENSOR: record entire field -----------------------
sensor.mask = true(Nx, Ny);
sensor.record = {'p'};

% -------------------- RUN SIMULATION ------------------------------------
input_args = {'DataCast', 'single', 'PMLSize', 20, 'PlotPML', false, 'PlotScale', [-0.5,0.5]};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
t_sim_elapsed=toc
%% 
% ---------- SETUP ----------
output_mp4 = true;    % set false to skip mp4
output_gif = false;    % set false to skip gif
mp4_filename = 'kwave_replay4.mp4';
gif_filename = 'kwave_replay3.gif';
fps = 25;             % frames per second for output
frame_delay = 1 / fps;

% prepare figure
hf = figure('Visible','off','Position',[100 100 600 600]);
colormap('parula');

% video writer
if output_mp4
    vw = VideoWriter(mp4_filename, 'MPEG-4');
    vw.FrameRate = fps;
    open(vw);
end

% compute color limits once for consistent coloring
clim = max(abs(sensor_data.p(:)));
if clim == 0, clim = 1; end

tic
% ---------- PREPARE COLOR LIMIT ----------
clim = max(abs(sensor_data.p(:)));
if clim == 0, clim = 1; end

% ---------- PRECOMPUTE FRAMES IN PARALLEL ----------
fprintf('Generating frames in parallel...\n');
frames = cell(1, Nt);

parfor t = 1:Nt
    pressure_snapshot = reshape(sensor_data.p(:,t), [Nx, Ny]);
    % Normalize to [0,1] for colormap
    normalized = (pressure_snapshot + clim) / (2*clim);
    normalized = min(max(normalized, 0), 1);
    % Apply parula colormap
    cmap = parula(256);
    idx = max(1, round(normalized * 255)); % 1–256
    rgb = ind2rgb(idx, cmap);
    frames{t} = im2uint8(rgb);
end

fprintf('Frame generation complete. Writing video...\n');

% ---------- WRITE MP4 SEQUENTIALLY ----------
if output_mp4
    vw = VideoWriter(mp4_filename, 'MPEG-4');
    vw.FrameRate = fps;
    open(vw);
    for t = 1:Nt
        writeVideo(vw, frames{t});
    end
    close(vw);
    fprintf('Saved MP4: %s\n', mp4_filename);
end

% ---------- WRITE GIF SEQUENTIALLY ----------
if output_gif
    for t = 1:Nt
        [A, map] = rgb2ind(frames{t}, 256);
        if t == 1
            imwrite(A, map, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', frame_delay);
        else
            imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
        end
    end
    fprintf('Saved GIF: %s\n', gif_filename);
end
t_render_elapsed=toc
