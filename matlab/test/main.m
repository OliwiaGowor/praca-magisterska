% === main.m: Main driver ===
clear; clc; close all;

% --- Path configuration ---
currentPath = fileparts(mfilename('fullpath'));
addpath(fullfile(currentPath, 'metrics'));
addpath(fullfile(currentPath, 'AES'));
addpath(fullfile(currentPath, 'DES'));
addpath(fullfile(currentPath, 'chaos_circular'));
addpath(fullfile(currentPath, 'blowfish'));
addpath(fullfile(currentPath, 'chacha'));
addpath(fullfile(currentPath, 'chaos_maps'));
addpath(fullfile(currentPath, 'chaos_bit_shuffling'));
addpath(fullfile(currentPath, 'chaos_2d'));
addpath(fullfile(currentPath, 'chaos_logistic'));
addpath(fullfile(currentPath, 'hyperchaotic_adaptive'));
addpath(fullfile(currentPath, 'images'));

% --- Test configuration ---
config.N_RUNS = 1;
run_attacks = true; % Set to false to skip
should_run_benchmarks = false; % Set to false to skip

% List of algorithm titles
config.titles = { ...
    'DES-CBC (C)', ...
    'AES-CBC (C)', ...
    'Blowfish-CBC (C)', ...
    'ChaCha20 (C)', ...
    'DES-CBC (MATLAB*)', ...
    'AES-CBC (MATLAB*)', ...
    'Blowfish-CBC (MATLAB*)', ...
    'ChaCha20 (MATLAB*)', ...
    'Chaos Circular Shifts (MATLAB)', ...
    '3D-Logistic-Chirikov (MATLAB)', ...
    'Chaos Bit Shuffling (MATLAB - Orig)', ...
    'Chaos Bit Shuffling (MATLAB* - Opt)', ...
    'Hyperchaotic and 2D Sensing (MATLAB*)', ...
    'Hu & Tian - Two-Stage Logistic (MATLAB*)', ...
    'Hyperchaotic Adaptive (MATLAB*)', ...
    'AES (Pixel/ECB Mode)', ...
    'DES (Pixel/ECB Mode)', ...
    };

fprintf('=== START TESTÓW (Przebiegi: %d) ===\n', config.N_RUNS);

[data] = setup_test_data();

if should_run_benchmarks
    raw_results = run_benchmarks(data, config);

    avg_results = process_results(raw_results, data, config);

    visualize_results(avg_results, raw_results, data, config);

end
% --- Attack Analysis ---

if run_attacks
    % Ensure metrics path is added
    addpath(fullfile(currentPath, 'metrics'));

    fprintf('\n--- Uruchamianie rozszerzonej analizy ataków (Alteration Attacks) ---\n');
    % Run comprehensive alteration analysis
    run_alteration_analysis(data);
end

fprintf('\n=== ZAKOŃCZONO POMYŚLNIE ===\n');
