function run_alteration_analysis(data)
    % RUN_ALTERATION_ANALYSIS Comprehensive Alteration Attack Analysis.
    %
    % Tests robustness of 15 algorithms (matching main.m order).
    % Saves results to CSV, MAT, and PNG files.
    
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    results_dir = fullfile(pwd, 'results_alteration');
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    
    fprintf('\n=== ROZPOCZYNAM KOMPLEKSOWĄ ANALIZĘ ATAKÓW (15 ALGORYTMÓW) ===\n');
    fprintf('Wyniki zostaną zapisane w: %s\n', results_dir);
    
    % --- Global Configuration ---
    target_size = data.size;
    target_numel = prod(target_size);
    
    % Generate keys for standard algorithms
    key_32 = uint8(randi([0 255], 1, 32)); 
    key_16 = uint8(randi([0 255], 1, 16));
    key_8  = uint8(randi([0 255], 1, 8));
    iv_16  = uint8(randi([0 255], 1, 16));
    iv_8   = uint8(randi([0 255], 1, 8));
    iv_12  = uint8(randi([0 255], 1, 12)); 
    key_des_hex = '1234567890ABCDEF'; 
    iv_des_hex  = 'FEDCBA0987654321';

    % --- Algorithm Definitions ---
    algos = {};
    
    % 1. DES-CBC (C)
    algos{end+1} = struct('name', 'DES-CBC (C)', ...
        'enc', @(img) crypto_mex(img(:)', key_8, iv_8, 'des-cbc', true), ...
        'dec', @(ct)  decrypt_and_restore(ct, @(c) crypto_mex(c, key_8, iv_8, 'des-cbc', false), target_size));

    % 2. AES-CBC (C)
    algos{end+1} = struct('name', 'AES-CBC (C)', ...
        'enc', @(img) crypto_mex(img(:)', key_32, iv_16, 'aes-256-cbc', true), ...
        'dec', @(ct)  decrypt_and_restore(ct, @(c) crypto_mex(c, key_32, iv_16, 'aes-256-cbc', false), target_size));

    % 3. Blowfish-CBC (C)
    algos{end+1} = struct('name', 'Blowfish-CBC (C)', ...
        'enc', @(img) crypto_mex(img(:)', key_16, iv_8, 'bf-cbc', true), ...
        'dec', @(ct)  decrypt_and_restore(ct, @(c) crypto_mex(c, key_16, iv_8, 'bf-cbc', false), target_size));

    % 4. ChaCha20 (C)
    algos{end+1} = struct('name', 'ChaCha20 (C)', ...
        'enc', @(img) crypto_mex(img(:)', key_32, iv_16, 'chacha20', true), ...
        'dec', @(ct)  decrypt_and_restore(ct, @(c) crypto_mex(c, key_32, iv_16, 'chacha20', false), target_size));

    % 5. DES-CBC (MATLAB)
    algos{end+1} = struct('name', 'DES-CBC (MATLAB)', ...
        'enc', @(img) des_enc_helper(img, key_des_hex, iv_des_hex), ...
        'dec', @(ct)  DES_ImgProcess_Opt.decrypt(ct, key_des_hex, iv_des_hex, target_size));

    % 6. AES-CBC (MATLAB)
    algos{end+1} = struct('name', 'AES-CBC (MATLAB)', ...
        'enc', @(img) aes_mode(img(:)', key_32, iv_16, 'encrypt', 'cbc'), ...
        'dec', @(ct)  decrypt_and_restore(ct, @(c) aes_mode(c, key_32, iv_16, 'decrypt', 'cbc'), target_size));

    % 7. Blowfish (MATLAB)
    bf_engine = ImageBlowfish(key_16);
    algos{end+1} = struct('name', 'Blowfish (MATLAB)', ...
        'enc', @(img) bf_engine.encryptImage(img, iv_8), ... 
        'dec', @(ct)  bf_engine.decryptImage(ct, iv_8, target_size, 'uint8'));

    % 8. ChaCha20 (MATLAB)
    algos{end+1} = struct('name', 'ChaCha20 (MATLAB)', ...
        'enc', @(img) chacha20(key_32, iv_12, 0, img(:)), ...
        'dec', @(ct)  reshape(chacha20(key_32, iv_12, 0, ct), target_size));

    % 9. Chaos Circular
    algos{end+1} = struct('name', 'Chaos Circular', 'type', 'chaos_stateful', ...
        'enc_gen', @(img) encrypt_circular_chaotic(img), ...
        'dec_gen', @(ct, keys) decrypt_circular_chaotic(ct, keys));

    % 10. 3D-Logistic-Chirikov
    algos{end+1} = struct('name', '3D-Logistic-Chirikov', 'type', 'chaos_stateful', ...
        'enc_gen', @(img) chaos_map_enc_helper(img), ...
        'dec_gen', @(ct, state) chaos_map_dec_helper(ct, state, target_size));

    % 11. Chaos Bit Shuffle (Orig)
    algos{end+1} = struct('name', 'Chaos Bit Shuffle (Orig)', 'type', 'chaos_stateful', ...
        'enc_gen', @(img) bit_shuffling_enc_wrapper(img, 'orig'), ... 
        'dec_gen', @(ct, keys) bit_shuffling_dec_wrapper(ct, keys, 'orig'));

    % 12. Chaos Bit Shuffle (Opt)
    algos{end+1} = struct('name', 'Chaos Bit Shuffle (Opt)', 'type', 'chaos_stateful', ...
        'enc_gen', @(img) bit_shuffling_enc_wrapper(img, 'opt'), ... 
        'dec_gen', @(ct, keys) bit_shuffling_dec_wrapper(ct, keys, 'opt'));

    % 13. Hyperchaotic 2D
    algos{end+1} = struct('name', 'Hyperchaotic 2D', 'type', 'chaos_stateful', ...
        'enc_gen', @(img) chaos_2d_adapter('encrypt', img), ...
        'dec_gen', @(ct, keys) chaos_2d_adapter('decrypt', ct, keys));

    % 14. Hu & Tian
    algos{end+1} = struct('name', 'Hu & Tian', 'type', 'chaos_stateful', ...
        'enc_gen', @(img) hu_tian_adapter('encrypt', img), ...
        'dec_gen', @(ct, keys) hu_tian_adapter('decrypt', ct, keys));

    % 15. Hyperchaotic Adaptive
    algos{end+1} = struct('name', 'Hyperchaotic Adaptive', 'type', 'chaos_stateful', ...
        'enc_gen', @(img) hyperchaotic_adaptive_wrapper(img, 'encrypt'), ...
        'dec_gen', @(ct, keys) hyperchaotic_adaptive_wrapper(ct, 'decrypt', keys));

    % --- 3. Attack Configuration ---
    attack_configs = {
        'contiguous_crop', [0.1, 0.25, 0.5, 0.75];
        'crop',            [0.1, 0.25, 0.5, 0.75];
        'tamper',          [0.1, 0.25, 0.5, 0.75];
        'salt_pepper',     [0.05, 0.10, 0.25, 0.5];
        'gaussian',        [0.00001, 0.0001, 0.001, 0.01];
    };
    
    % --- Load Tamper Image (Pattern) ---
    try
        img_tamper = imread('peppers.png');
        if size(img_tamper,3)>1, img_tamper=rgb2gray(img_tamper); end
        img_tamper = imresize(img_tamper, target_size);
    catch
        % Fallback: Generate a Checkerboard Pattern instead of Random Noise
        % (More distinct for tampering)
        [X, Y] = meshgrid(1:target_size(2), 1:target_size(1));
        img_tamper = uint8((mod(floor(X/20) + floor(Y/20), 2) == 0) * 255);
    end

    % --- 4. Main Testing Loop ---
    results = struct();
    visual_data = cell(length(algos), size(attack_configs, 1));

    for a = 1:length(algos)
        algo = algos{a};
        fprintf('[%d/%d] Algorytm: %-30s ... ', a, length(algos), algo.name);
        
        try
            if isfield(algo, 'type') && strcmp(algo.type, 'chaos_stateful')
                [C_base, state] = algo.enc_gen(data.img_orig);
                decrypt_fn = @(ct) algo.dec_gen(ct, state);
            else
                C_base = algo.enc(data.img_orig);
                decrypt_fn = algo.dec;
            end
        catch ME
            fprintf('NIEPOWODZENIE (Szyfr): %s\n', ME.message);
            continue;
        end
        
        for t = 1:size(attack_configs, 1)
            att_type = attack_configs{t, 1};
            intensities = attack_configs{t, 2};
            
            psnr_vec = zeros(1, length(intensities));
            ssim_vec = zeros(1, length(intensities));
            sd_vec   = zeros(1, length(intensities));
            
            vis_results = struct('C_img', {}, 'P_img', {}, 'intensity', {});
            
            for i = 1:length(intensities)
                inte = intensities(i);
                
                % --- MODYFIKACJA: Ataki w LEWYM GÓRNYM ROGU (1,1) ---
                if strcmp(att_type, 'tamper') || strcmp(att_type, 'contiguous_crop')
                    C_prim = C_base;
                    
                    % 1. Oblicz ile pikseli zmienić
                    total_pixels = numel(C_base);
                    affected_pixels = round(total_pixels * inte);
                    
                    if isvector(C_prim)
                        % Wektor: Zmieniamy początek (indeks 1)
                        end_idx = min(total_pixels, affected_pixels);
                        
                        if strcmp(att_type, 'tamper')
                            img_t_flat = img_tamper(:);
                            % Zabezpieczenie na wypadek gdyby img_tamper był mniejszy
                            len_to_copy = min(numel(img_t_flat), end_idx);
                            C_prim(1:len_to_copy) = img_t_flat(1:len_to_copy);
                        else
                            % Contiguous Crop -> Wyzeruj początek
                            C_prim(1:end_idx) = 0;
                        end
                    else
                        % Macierz 2D: Kwadrat w lewym górnym rogu
                        block_side = round(sqrt(affected_pixels));
                        [rows, cols] = size(C_prim);
                        
                        % Ustal granice bloku w rogu (1,1)
                        r_end = min(rows, block_side);
                        c_end = min(cols, block_side);
                        
                        if strcmp(att_type, 'tamper')
                            C_prim(1:r_end, 1:c_end) = img_tamper(1:r_end, 1:c_end);
                        else
                            % Contiguous Crop -> Wyzeruj blok w rogu
                            C_prim(1:r_end, 1:c_end) = 0;
                        end
                    end
                else
                    % Pozostałe ataki (szum, losowe wycinanie) obsługiwane standardowo
                    C_prim = apply_attack(C_base, att_type, inte, img_tamper);
                end
                
                % Wyciszanie ostrzeżeń
                w_state = warning('off', 'all'); 
                try
                    P_prim = decrypt_fn(C_prim);
                catch
                    P_prim = zeros(target_size, 'uint8');
                end
                warning(w_state);
                
                if ~isequal(size(P_prim), target_size)
                    P_tmp = zeros(target_size, 'uint8');
                    num_copy = min(numel(P_prim), target_numel);
                    P_prim_flat = P_prim(:);
                    P_tmp(1:num_copy) = P_prim_flat(1:num_copy);
                    P_prim = P_tmp;
                end
                
                [p, s, d] = calculate_alteration_metrics(data.img_orig, P_prim);
                psnr_vec(i) = p; ssim_vec(i) = s; sd_vec(i) = d;
                
                if isvector(C_prim)
                     C_disp = zeros(target_size, 'uint8');
                     n_c = min(numel(C_prim), target_numel);
                     C_disp(1:n_c) = C_prim(1:n_c);
                else
                     C_disp = C_prim;
                end
                vis_results(i).C_img = C_disp;
                vis_results(i).P_img = P_prim;
                vis_results(i).intensity = inte;
            end
            visual_data{a, t} = vis_results;
            
            key_res = sprintf('a%d_%s', a, att_type);
            results.(key_res).psnr = psnr_vec;
            results.(key_res).ssim = ssim_vec;
            results.(key_res).sd   = sd_vec;
            results.(key_res).intensities = intensities;
        end
        fprintf('Gotowe.\n');
    end
    
    % --- 5. Saving Results ---
    fprintf('Zapisywanie wyników do plików...\n');
    mat_filename = fullfile(results_dir, ['alteration_results_' timestamp '.mat']);
    save(mat_filename, 'results', 'algos', 'attack_configs', 'visual_data');
    fprintf('  -> Dane MAT: %s\n', mat_filename);

    csv_filename = fullfile(results_dir, ['alteration_metrics_' timestamp '.csv']);
    fid = fopen(csv_filename, 'w');
    if fid == -1
        warning('Nie można otworzyć pliku CSV do zapisu.');
    else
        fprintf(fid, 'Algorithm,Attack_Type,Intensity,PSNR,SSIM,SD\n');
        for a = 1:length(algos)
            algo_name = algos{a}.name;
            for t = 1:size(attack_configs, 1)
                att_type = attack_configs{t, 1};
                key_res = sprintf('a%d_%s', a, att_type);
                if isfield(results, key_res)
                    res = results.(key_res);
                    ints = res.intensities;
                    for k = 1:length(ints)
                        fprintf(fid, '%s,%s,%.6f,%.6f,%.6f,%.6f\n', ...
                            algo_name, att_type, ints(k), ...
                            res.psnr(k), res.ssim(k), res.sd(k));
                    end
                end
            end
        end
        fclose(fid);
        fprintf('  -> Dane CSV: %s\n', csv_filename);
    end
    
    % --- 6. Visualization & Plot Saving ---
    fprintf('Generowanie i zapisywanie wykresów...\n');
    plot_alteration_results(results, algos, attack_configs, results_dir, timestamp);
    plot_attack_examples(visual_data, algos, attack_configs, results_dir, timestamp);
    
    fprintf('=== ANALIZA ZAKOŃCZONA ===\n');
end

% --- Helper Functions ---
function restore_warnings(w1, w2)
    warning(w1);
    warning(w2);
end

function img_out = decrypt_and_restore(ct_flat, dec_core_fn, target_size)
    pt_padded = dec_core_fn(ct_flat);
    target_numel = prod(target_size);
    if numel(pt_padded) >= target_numel
        pt_truncated = pt_padded(1:target_numel);
    else
        pt_truncated = zeros(1, target_numel, 'uint8');
        pt_truncated(1:numel(pt_padded)) = pt_padded;
    end
    img_out = reshape(pt_truncated, target_size);
end

function out = des_enc_helper(img, k, v)
    [ct, ~] = DES_ImgProcess_Opt.encrypt(img, k, v);
    out = ct;
end

function [C, keys] = bit_shuffling_enc_wrapper(img, mode)
    if strcmp(mode, 'opt')
        [C_double, keys] = chaos_bit_shuffling_optimized('encrypt', img);
    else
        if exist('chaos_bit_shuffling', 'file') == 2
             [C_double, keys] = chaos_bit_shuffling('encrypt', img);
        else
             warning('Original bit shuffling not found, using optimized.');
             [C_double, keys] = chaos_bit_shuffling_optimized('encrypt', img);
        end
    end
    C = uint8(C_double);
end

function img_out = bit_shuffling_dec_wrapper(ct, keys, mode)
    if strcmp(mode, 'opt')
        img_out = chaos_bit_shuffling_optimized('decrypt', ct, keys);
    else
        if exist('chaos_bit_shuffling', 'file') == 2
             img_out = chaos_bit_shuffling('decrypt', ct, keys);
        else
             img_out = chaos_bit_shuffling_optimized('decrypt', ct, keys);
        end
    end
    img_out = uint8(img_out);
end

function [out, state] = hyperchaotic_adaptive_wrapper(img, mode, state)
    if exist('hc_key_gen', 'file') ~= 2
        if exist('hyperchaotic_adaptive', 'dir')
            addpath('hyperchaotic_adaptive');
        elseif exist(fullfile('test1', 'hyperchaotic_adaptive'), 'dir')
            addpath(fullfile('test1', 'hyperchaotic_adaptive'));
        end
    end
    if strcmp(mode, 'encrypt')
        if exist('hc_key_gen', 'file') == 2 && exist('hc_encrypt', 'file') == 2
            [seq, params] = hc_key_gen(img); 
            out = hc_encrypt(img, seq);
            state.seq = seq;
            state.params = params; 
        else
            warning('HC Adaptive functions not found. Using dummy.');
            out = img; 
            state.seq = [];
        end
    else
        if exist('hc_decrypt', 'file') == 2 && ~isempty(state) && isfield(state, 'seq')
            out = hc_decrypt(img, state.seq, size(img));
        else
             out = img;
        end
    end
end

function [ct, state] = chaos_map_enc_helper(img)
    [rows, cols, ~] = size(img);
    num_pixels = rows * cols;
    key_chirikov = rand(1, 4) * 10;
    [diffused_img, seqs] = logistic_diffusion_core(img, num_pixels);
    [ct, ~] = encrypt(diffused_img, key_chirikov);
    state.key = key_chirikov;
    state.seqs = seqs;
end

function img_out = chaos_map_dec_helper(ct, state, sz)
    dec_chirikov = decrypt(ct, state.key);
    img_out = logistic_diffusion_inv_core(dec_chirikov, state.seqs, sz);
end

function [out_img, seqs] = logistic_diffusion_core(img, N)
    x = zeros(1, N+1); y = zeros(1, N+1); z = zeros(1, N+1);
    x(1)=0.2350; y(1)=0.3500; z(1)=0.7350;
    a=0.0125; b=0.0157; l=3.7700;
    for i=1:N
        x(i+1) = l*x(i)*(1-x(i)) + b*y(i)*y(i)*x(i) + a*z(i)*z(i)*z(i);
        y(i+1) = l*y(i)*(1-y(i)) + b*z(i)*z(i)*y(i) + a*x(i)*x(i)*x(i);
        z(i+1) = l*z(i)*(1-z(i)) + b*x(i)*x(i)*z(i) + a*y(i)*y(i)*y(i);
    end
    x = x(2:end); y = y(2:end); z = z(2:end);
    seqs.Sx = uint8(ceil(mod((x*1000000), 256)));
    seqs.Sy = uint8(ceil(mod((y*1000000), 256)));
    seqs.Sz = uint8(ceil(mod((z*1000000), 256)));
    seqs.x=x; seqs.y=y; seqs.z=z;
    [r, c, ch] = size(img);
    PR = double(reshape(img(:,:,1), 1, []));
    CDR = x .* PR; 
    CCR = bitxor(seqs.Sx, uint8(CDR));
    out_img = reshape(CCR, r, c);
end

function out_img = logistic_diffusion_inv_core(img, seqs, sz)
    [r, c, ch] = size(img);
    xi = 1./seqs.x;
    DR = reshape(img(:,:,1), 1, []);
    DDR = bitxor(seqs.Sx, uint8(DR));
    DDDR = xi .* double(DDR);
    out_img = uint8(reshape(DDDR, r, c));
end

% --- Visualization Functions ---
function plot_alteration_results(results, algos, attack_configs, save_dir, timestamp)
    metrics = {'psnr', 'ssim', 'sd'};
    metric_titles = {'PSNR (Peak Signal-to-Noise Ratio)', 'SSIM (Structural Similarity)', 'SD (Spectral Distortion)'};
    num_attacks = size(attack_configs, 1);
    num_metrics = 3;

    % --- PALETA (15 kolorów) ---
    colors = [
        0.0000, 0.0000, 0.0000; % #000000
        0.0000, 0.1176, 0.6392; % #001ea3
        0.2314, 0.3059, 0.8902; % #3b4ee3
        0.4392, 0.5059, 1.0000; % #7081ff
        0.4706, 0.2941, 0.8235; % #784bd2
        0.6196, 0.3882, 0.8118; % #9e63cf
        0.6706, 0.2392, 0.7216; % #ab3db8
        0.6588, 0.1412, 0.2784; % #a82447
        0.8784, 0.2431, 0.4667; % #e03e77
        0.9608, 0.4784, 0.5373; % #f57a89
        1.0000, 0.5176, 0.0000; % #ff8400
        1.0000, 0.7020, 0.1020; % #ffb31a
        0.9490, 1.0000, 0.0000; % #f2ff00
        0.6157, 0.7098, 0.0000; % #9db500
        0.4157, 0.4784, 0.0000  % #6a7a00
    ];

    linestyles = {'-', '--', ':', '-.'};
    if length(algos) > size(colors, 1)
        colors = [colors; rand(length(algos)-size(colors,1), 3)];
    end
    
    for t = 1:num_attacks
        att_type = attack_configs{t, 1};
        switch att_type
            case 'contiguous_crop', pretty_att = 'Wycięcie Ciągłe'; safe_att = 'ContiguousCrop';
            case 'crop',            pretty_att = 'Wycięcie Losowe'; safe_att = 'RandomCrop';
            case 'tamper',          pretty_att = 'Podmiana (Tamper)'; safe_att = 'Tamper';
            case 'salt_pepper',     pretty_att = 'Szum Sól/Pieprz'; safe_att = 'SaltPepper';
            case 'gaussian',        pretty_att = 'Szum Gaussa'; safe_att = 'Gaussian';
            otherwise,              pretty_att = att_type; safe_att = att_type;
        end
        
        for m = 1:num_metrics
            metric_key = metrics{m};
            pretty_metric = metric_titles{m};
            f = figure('Name', [pretty_att ' - ' metric_key], 'Color', 'w', 'Position', [100 100 1000 700], 'Visible', 'off');
            hold on;
            plot_handles = [];
            for a = 1:length(algos)
                key_res = sprintf('a%d_%s', a, att_type);
                if isfield(results, key_res)
                    res = results.(key_res);
                    color_idx = mod(a-1, size(colors,1)) + 1;
                    style_idx = mod(a-1, length(linestyles)) + 1;
                    fmt = [linestyles{style_idx} 'o'];
                    p = plot(res.intensities, res.(metric_key), fmt, ...
                         'Color', colors(color_idx, :), 'LineWidth', 2.0, 'MarkerSize', 5, 'DisplayName', algos{a}.name);
                    plot_handles = [plot_handles, p];
                end
            end
            hold off; grid on;
            title(sprintf('%s | %s\nData: %s', pretty_att, pretty_metric), 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'none');
            xlabel('Intensywność ataku', 'FontSize', 11);
            ylabel(['Wartość ' metric_key], 'FontSize', 11);
            if ~isempty(plot_handles)
                lgd = legend(plot_handles, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'none');
                try set(lgd, 'NumColumns', 4); catch; end
            end
            fname = fullfile(save_dir, sprintf('%s_%s_%s.png', metric_key, safe_att, timestamp));
            print(f, fname, '-dpng', '-r150'); 
            close(f);
        end
    end
end

function plot_attack_examples(visual_data, algos, attack_configs, save_dir, timestamp)
    num_attacks = size(attack_configs, 1);
    for a = 1:length(algos)
        for t = 1:num_attacks
            vis_results = visual_data{a, t};
            if isempty(vis_results), continue; end
            num_intensities = length(vis_results);
            att_type = attack_configs{t, 1};
            switch att_type
                case 'contiguous_crop', pretty = 'Wycięcie Ciągłe'; safe_att = 'ContiguousCrop';
                case 'crop',            pretty = 'Wycięcie Losowe'; safe_att = 'RandomCrop';
                case 'tamper',          pretty = 'Podmiana (Tamper)'; safe_att = 'Tamper';
                case 'salt_pepper',     pretty = 'Szum Sól/Pieprz'; safe_att = 'SaltPepper';
                case 'gaussian',        pretty = 'Szum Gaussa'; safe_att = 'Gaussian';
                otherwise,              pretty = att_type; safe_att = att_type;
            end
            safe_algo = regexprep(algos{a}.name, '[^a-zA-Z0-9]', '_');
            f = figure('Name', [algos{a}.name ' - ' pretty], 'Color', 'w', ...
                'Position', [100, 100, 300 * num_intensities, 600], 'Visible', 'off');
            for k = 1:num_intensities
                data_struct = vis_results(k);
                subplot(2, num_intensities, k);
                imshow(data_struct.C_img);
                title(sprintf('Intens.: %.2f\nSzyfrogram', data_struct.intensity), 'FontSize', 10);
                subplot(2, num_intensities, k + num_intensities);
                imshow(data_struct.P_img);
                title('Odtworzony', 'FontSize', 10);
            end
            sgtitle(sprintf('%s: %s', algos{a}.name, pretty), 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'none');
            fname = fullfile(save_dir, sprintf('Vis_%s_%s_%s.png', safe_att, safe_algo, timestamp));
            print(f, fname, '-dpng', '-r150');
            close(f);
        end
    end
    fprintf('  -> Zapisano wizualizacje ataków (wszystkie intensywności).\n');
end