function res = run_benchmarks(data, config)
    % RUN_BENCHMARKS - Wykonuje testy (W tym 3-punktowy test wrażliwości)
    
    N = config.N_RUNS;
    num_algos = length(config.titles);
    
    res = init_results_struct(num_algos, N);
    
    % Lista funkcji wrapperów dla algorytmów (kolejność musi pasować do config.titles)
    % Uwaga: Musisz upewnić się, że kolejność tutaj odpowiada main.m!
    algos = {
        @algo_des_c, ...
        @algo_aes_c, ...
        @algo_blowfish_c, ...
        @algo_chacha_c, ...
        @algo_des_matlab, ...
        @algo_aes_matlab, ...
        @algo_blowfish_matlab, ...
        @algo_chacha_matlab, ...
        @algo_chaos_circular_matlab, ...
        @algo_chaos_map_matlab, ...
        @algo_chaos_bit_shuffling_original, ...
        @algo_chaos_bit_shuffling, ... % Optimized
        @algo_hyperchaotic_2d_matlab, ...
        @algo_hu_tian, ...
        @algo_hyperchaotic_adaptive, ...
        @algo_aes_pixel, ...  % <-- DODAJ TU
        @algo_des_pixel, ...  % <-- DODAJ TU
    };
    
    for run_idx = 1:N
        fprintf('--- Przebieg %d / %d ---\n', run_idx, N);
        
        for i = 1:length(algos)
            algo_func = algos{i};
            algo_name = config.titles{i};
            
            try
                % --- 1. Test "Start" (Pixel 1) ---
                d_start = data; 
                d_start.img_mod = data.img_mod_start; 
                d_start.img_flat_mod = data.img_flat_mod_start;
                [t1, ~, C1_start, C2_start, ~] = algo_func(d_start);
                
                % --- 2. Test "Mid" (Pixel Center) ---
                d_mid = data; 
                d_mid.img_mod = data.img_mod_mid; 
                d_mid.img_flat_mod = data.img_flat_mod_mid;
                [t2, t_dec, C1_mid, C2_mid, PT] = algo_func(d_mid); % Tu pobieramy też t_dec i PT
                
                % --- 3. Test "End" (Pixel End) ---
                d_end = data; 
                d_end.img_mod = data.img_mod_end; 
                d_end.img_flat_mod = data.img_flat_mod_end;
                [t3, ~, C1_end, C2_end, ~] = algo_func(d_end);
                
                % --- Obliczanie Metryk ---
                
                % Czas (Średnia z 3 szyfrowań dla stabilności)
                res.times_enc(i, run_idx) = mean([t1, t2, t3]);
                res.times_dec(i, run_idx) = t_dec;
                
                % Jakość (Entropy/Corr liczymy dla środkowego wariantu)
                res.entropy(i, run_idx) = calculate_entropy(C1_mid);
                [ch, cv, cd] = calculate_correlation(C1_mid);
                res.corr_h(i, run_idx) = ch;
                res.corr_v(i, run_idx) = cv;
                res.corr_d(i, run_idx) = cd;
                
                % NPCR / UACI (Dla każdego wariantu osobno)
                [n_s, u_s] = calculate_npcr_uaci(C1_start, C2_start);
                [n_m, u_m] = calculate_npcr_uaci(C1_mid, C2_mid);
                [n_e, u_e] = calculate_npcr_uaci(C1_end, C2_end);
                
                res.npcr_start(i, run_idx) = n_s; res.uaci_start(i, run_idx) = u_s;
                res.npcr_mid(i, run_idx)   = n_m; res.uaci_mid(i, run_idx)   = u_m;
                res.npcr_end(i, run_idx)   = n_e; res.uaci_end(i, run_idx)   = u_e;
                
                % Zapisz obrazy (dla galerii - wariant Mid)
                if run_idx == 1
                    res.images_enc{i} = C1_mid;
                    res.images_dec{i} = PT;
                end
                
            catch ME
                warning('%s error: %s', algo_name, ME.message);
                % fprintf('%s\n', getReport(ME));
            end
        end
    end
end

function res = init_results_struct(rows, cols)
    res.times_enc = NaN(rows, cols);
    res.times_dec = NaN(rows, cols);
    
    % Standard metrics
    res.entropy = NaN(rows, cols);
    res.corr_h = NaN(rows, cols);
    res.corr_v = NaN(rows, cols);
    res.corr_d = NaN(rows, cols);
    
    % NEW: Separate NPCR/UACI for 3 positions
    res.npcr_start = NaN(rows, cols); res.uaci_start = NaN(rows, cols);
    res.npcr_mid   = NaN(rows, cols); res.uaci_mid   = NaN(rows, cols);
    res.npcr_end   = NaN(rows, cols); res.uaci_end   = NaN(rows, cols);
    
    res.images_enc = cell(rows, 1);
    res.images_dec = cell(rows, 1);
end