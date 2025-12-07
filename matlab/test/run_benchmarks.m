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
