function res = run_benchmarks(data, config)
    % Main driver - executes individual algorithms
    
    N = config.N_RUNS;
    num_algos = length(config.titles);
    
    res = init_results_struct(num_algos, N);
    
    for run_idx = 1:N
        fprintf('--- Przebieg %d / %d ---\n', run_idx, N);
        
        % --- 1. AES (C) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_aes_c(data);
            res = helper_update_metrics(res, 1, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME, warning('AES (C) error: %s', ME.message); end

        % --- 2. Blowfish (C) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_blowfish_c(data);
            res = helper_update_metrics(res, 2, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME, warning('Blowfish (C) error: %s', ME.message); end

        % --- 3. ChaCha20 (C) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_chacha_c(data);
            res = helper_update_metrics(res, 3, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME, warning('ChaCha (C) error: %s', ME.message); end

        % --- 4. Chaos (MATLAB) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_chaos_circular_matlab(data);
            res = helper_update_metrics(res, 4, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME, warning('Chaos (MATLAB) error: %s', ME.message); end

        % --- 5. DES (MATLAB) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_des_matlab(data);
            res = helper_update_metrics(res, 5, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME, warning('DES (MATLAB) error: %s', ME.message); end

        % --- 6. AES (MATLAB) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_aes_matlab(data);
            res = helper_update_metrics(res, 6, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME, warning('AES (MATLAB) error: %s', ME.message); end

        % --- 7. Blowfish (MATLAB) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_blowfish_matlab(data);
            res = helper_update_metrics(res, 7, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME
            warning('Blowfish (MATLAB) error: %s', ME.message); 
            % Print details because Blowfish.m is sensitive to missing constants file
            fprintf('%s\n', getReport(ME));
        end
        
        % --- 8. ChaCha20 (Pure MATLAB) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_chacha_matlab(data);
            res = helper_update_metrics(res, 8, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME
            warning('ChaCha20 (MATLAB) error: %s', ME.message); 
            fprintf('%s\n', getReport(ME));
        end
        
        % --- 9. 3D Logistic + Chirikov (MATLAB) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_chaos_map_matlab(data);
            res = helper_update_metrics(res, 9, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME
            warning('Chaotic Map (MATLAB) error: %s', ME.message);
            fprintf('%s\n', getReport(ME));
        end
        
        % --- 10. Chaos Bit Shuffling (Moysis) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_chaos_bit_shuffling(data);
            res = helper_update_metrics(res, 10, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME
            warning('Chaos Bit Shuffling error: %s', ME.message);
            % Opcjonalnie wypisz pełny błąd dla debugowania, bo ten algorytm jest skomplikowany
            % fprintf('%s\n', getReport(ME));
        end
        
        % % --- 11. Hyperchaotic 2D (S. Liu et al.) ---
        try
            % index 11 musi odpowiadać pozycji w config.titles
            [t_enc, t_dec, C1, C2, PT] = algo_hyperchaotic_2d_matlab(data);
            res = helper_update_metrics(res, 11, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME
            warning('Hyperchaotic 2D error: %s', ME.message);
            fprintf('%s\n', getReport(ME));
        end
        
        % --- 12. Hu & Tian (Two-Stage Logistic) ---
        try
            [t_enc, t_dec, C1, C2, PT] = algo_hu_tian(data);
            res = helper_update_metrics(res, 12, run_idx, t_enc, t_dec, C1, C2, PT);
        catch ME
            warning('Hu & Tian error: %s', ME.message);
            % fprintf('%s\n', getReport(ME)); 
        end
        
    end
end

function res = init_results_struct(rows, cols)
    res.times_enc = NaN(rows, cols);
    res.times_dec = NaN(rows, cols);
    res.npcr = NaN(rows, cols);
    res.uaci = NaN(rows, cols);
    res.entropy = NaN(rows, cols);
    res.corr_h = NaN(rows, cols);
    res.corr_v = NaN(rows, cols);
    res.corr_d = NaN(rows, cols);
    res.images_enc = cell(rows, 1);
    res.images_dec = cell(rows, 1);
end
