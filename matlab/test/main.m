% === main.m: Main driver (NO WARM-UP) ===
clear; clc; close all;

% --- Path configuration ---
currentPath = fileparts(mfilename('fullpath'));
addpath(fullfile(currentPath, 'metrics'));
addpath(fullfile(currentPath, 'AES')); 
addpath(fullfile(currentPath, 'DES'));
addpath(fullfile(currentPath, 'chaos_circular'));
addpath(fullfile(currentPath, 'blowfish'));
addpath(fullfile(currentPath, 'chacha'));

% --- Test configuration ---
config.N_RUNS = 2; 

% List of algorithm titles
config.titles = { ...
    'AES-CBC (Hardware-C)', ...
    'Blowfish-CBC (Software-C)', ...
    'ChaCha20 (Software-C)', ...
    'Chaotyczny (MATLAB)', ...
    'DES-CBC (MATLAB)', ...
    'AES-CBC (Pure MATLAB)', ...
    'Blowfish-CBC (Pure MATLAB)', ...
    'ChaCha20 (Pure MATLAB)' ... % New Entry (Index 8)
};

fprintf('=== START TESTÓW (Przebiegi: %d) ===\n', config.N_RUNS);

% 1. Data preparation
[data] = setup_test_data();

% 2. Run benchmarks
raw_results = run_benchmarks(data, config);

% 3. Process and save results
avg_results = process_results(raw_results, data, config);

% 4. Visualization
visualize_results(avg_results, raw_results, data, config);

fprintf('\n=== ZAKOŃCZONO POMYŚLNIE ===\n');
