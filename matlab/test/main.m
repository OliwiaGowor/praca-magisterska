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
addpath(fullfile(currentPath, 'chaos_maps'));
addpath(fullfile(currentPath, 'chaos_bit_shuffling'));
addpath(fullfile(currentPath, 'chaos_2d'));
addpath(fullfile(currentPath, 'chaos_logistic'));

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
    'ChaCha20 (Pure MATLAB)', ...
    '3D-Logistic-Chirikov (MATLAB)', ...
    'Chaos Bit Shuffling (MATLAB)', ...
    'Hyperchaotic 2D (MATLAB*)', ...
    'Hu & Tian - Two-Stage Logistic (MATLAB*)', ...
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

% --- Analiza Ataków ---
run_attacks = true; % Ustaw na false jeśli chcesz pominąć (trwa dłużej)
if run_attacks
    % Upewnij się, że ścieżki są dodane
    addpath(fullfile(currentPath, 'metrics')); 
    
    fprintf('\n--- Uruchamianie rozszerzonej analizy ataków (Alteration Attacks) ---\n');
    % Wywołujemy nową funkcję przekazując dane (obrazek)
    run_alteration_analysis(data);
end

fprintf('\n=== ZAKOŃCZONO POMYŚLNIE ===\n');
