function startup()
% STARTUP Add the project folders and bundled k-Wave toolbox to the MATLAB path.

repo_root = fileparts(mfilename('fullpath'));

addpath(fullfile(repo_root, 'src'));
addpath(fullfile(repo_root, 'analysis'));
addpath(fullfile(repo_root, 'legacy'));
addpath(fullfile(repo_root, 'notebooks'));
addpath(genpath(fullfile(repo_root, 'k-wave-toolbox-version-1.4', 'k-Wave')));

fprintf('Project paths added from %s\n', repo_root);
end
