sizex_vals = [0.10, 0.20, 0.30, 0.40];
dx_vals    = [1e-3, (1e-3)/2, 1e-4];
f0_vals    = [2e4, 4e4, 6e4, 1e5];

[results, stats] = convergence_analysis2(sizex_vals, dx_vals, f0_vals);


%% 

sizex = 0.20;
dx    = (1e-3)/2;
f0    = 1e5;
simulate_plane('sizex', sizex, 'dx', dx, 'f0', f0);
