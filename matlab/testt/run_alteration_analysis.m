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
    
    % --- Suppress known warnings during attacks ---
    orig_warn_state_aes = warning('off', 'AES:BadPadding');
    orig_warn_state_des = warning('off', 'MATLAB:subsassigndimmismatch'); 
    cleanupObj = onCleanup(@() restore_warnings(orig_warn_state_aes, orig_warn_state_des)); 

    % --- 1. Global Configuration ---
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

    % --- 2. Algorithm Definitions ---
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
        'enc', @(img) crypto_mex(img(:)', key_32, iv_12, 'chacha20', true), ...
        'dec', @(ct)  decrypt_and_restore(ct, @(c) crypto_mex(c, key_32, iv_12, 'chacha20', false), target_size));

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

    % 15. Hyperchaotic Adaptive (POPRAWIONA IMPLEMENTACJA)
    algos{end+1} = struct('name', 'Hyperchaotic Adaptive', 'type', 'chaos_stateful', ...
        'enc_gen', @(img) hyperchaotic_adaptive_wrapper(img, 'encrypt'), ...
        'dec_gen', @(ct, keys) hyperchaotic_adaptive_wrapper(ct, 'decrypt', keys));

    % --- 3. Attack Configuration ---
    attack_configs = {
        'contiguous_crop', [0.1, 0.2, 0.3, 0.5];
        'crop',            [0.1, 0.2, 0.3, 0.5];
        'tamper',          [0.05, 0.1, 0.2, 0.4];
        'salt_pepper',     [0.01, 0.05, 0.1, 0.2];
        'gaussian',        [0.005, 0.01, 0.02, 0.05];
    };
    
    vis_intensities = containers.Map;
    vis_intensities('contiguous_crop') = 0.2;
    vis_intensities('crop')            = 0.2;
    vis_intensities('tamper')          = 0.1;
    vis_intensities('salt_pepper')     = 0.05;
    vis_intensities('gaussian')        = 0.01;

    try
        img_tamper = imread('peppers.png');
        if size(img_tamper,3)>1, img_tamper=rgb2gray(img_tamper); end
        img_tamper = imresize(img_tamper, target_size);
    catch
        img_tamper = uint8(randi([0 255], target_size));
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
            vis_int = 0.1; 
            if isKey(vis_intensities, att_type), vis_int = vis_intensities(att_type); end
            
            psnr_vec = zeros(1, length(intensities));
            ssim_vec = zeros(1, length(intensities));
            sd_vec   = zeros(1, length(intensities));
            
            for i = 1:length(intensities)
                inte = intensities(i);
                C_prim = apply_attack(C_base, att_type, inte, img_tamper);
                
                try
                    P_prim = decrypt_fn(C_prim);
                catch
                    P_prim = zeros(target_size, 'uint8');
                end
                
                if ~isequal(size(P_prim), target_size)
                    P_tmp = zeros(target_size, 'uint8');
                    num_copy = min(numel(P_prim), target_numel);
                    P_prim_flat = P_prim(:);
                    P_tmp(1:num_copy) = P_prim_flat(1:num_copy);
                    P_prim = P_tmp;
                end
                
                [p, s, d] = calculate_alteration_metrics(data.img_orig, P_prim);
                psnr_vec(i) = p; ssim_vec(i) = s; sd_vec(i) = d;
                
                if abs(inte - vis_int) < 0.001 || i == 1
                     if isvector(C_prim)
                         C_disp = zeros(target_size, 'uint8');
                         n_c = min(numel(C_prim), target_numel);
                         C_disp(1:n_c) = C_prim(1:n_c);
                     else
                         C_disp = C_prim;
                     end
                     vis_struct.C_img = C_disp;
                     vis_struct.P_img = P_prim;
                     vis_struct.intensity = inte;
                     visual_data{a, t} = vis_struct;
                end
            end
            
            key_res = sprintf('a%d_%s', a, att_type);
            results.(key_res).psnr = psnr_vec;
            results.(key_res).ssim = ssim_vec;
            results.(key_res).sd   = sd_vec;
            results.(key_res).intensities = intensities;
        end
        fprintf('Gotowe.\n');
    end
    
    % --- 5. Saving Results (CSV & MAT) ---
    fprintf('Zapisywanie wyników do plików...\n');
    
    % 5a. Save MAT file (Complete Workspace)
    mat_filename = fullfile(results_dir, ['alteration_results_' timestamp '.mat']);
    save(mat_filename, 'results', 'algos', 'attack_configs', 'visual_data');
    fprintf('  -> Dane MAT: %s\n', mat_filename);

    % 5b. Save CSV file (Metrics for Plotting)
    csv_filename = fullfile(results_dir, ['alteration_metrics_' timestamp '.csv']);
    fid = fopen(csv_filename, 'w');
    if fid == -1
        warning('Nie można otworzyć pliku CSV do zapisu.');
    else
        % Header
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

% --- UPDATED: Hyperchaotic Adaptive Wrapper ---
function [out, state] = hyperchaotic_adaptive_wrapper(img, mode, state)
    % 1. Ensure path exists
    if exist('hc_key_gen', 'file') ~= 2
        % Try standard paths
        if exist('hyperchaotic_adaptive', 'dir')
            addpath('hyperchaotic_adaptive');
        elseif exist(fullfile('test1', 'hyperchaotic_adaptive'), 'dir')
            addpath(fullfile('test1', 'hyperchaotic_adaptive'));
        end
    end

    if strcmp(mode, 'encrypt')
        % Encrypt
        if exist('hc_key_gen', 'file') == 2 && exist('hc_encrypt', 'file') == 2
            % Generate key/sequence based on image hash
            [seq, params] = hc_key_gen(img); 
            out = hc_encrypt(img, seq);
            
            % Save state for decryption
            state.seq = seq;
            state.params = params; 
        else
            warning('HC Adaptive functions not found (hc_key_gen/hc_encrypt). Using dummy.');
            out = img; 
            state.seq = [];
        end
    else
        % Decrypt
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
    f = figure('Name', 'Analiza Metryk', 'Color', 'w', 'Position', [50 50 1400 900]);
    metrics = {'psnr', 'ssim', 'sd'};
    titles = {'PSNR (Im wyżej tym lepiej)', 'SSIM (Im bliżej 1 tym lepiej)', 'SD (Im niżej tym lepiej)'};
    rows = 3; cols = size(attack_configs, 1);
    
    colors = turbo(length(algos)); 
    
    for m = 1:rows
        metric_key = metrics{m};
        for t = 1:cols
            att_type = attack_configs{t, 1};
            subplot(rows, cols, (m-1)*cols + t);
            hold on;
            for a = 1:length(algos)
                key_res = sprintf('a%d_%s', a, att_type);
                if isfield(results, key_res)
                    res = results.(key_res);
                    plot(res.intensities, res.(metric_key), '-o', ...
                         'Color', colors(a,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
                         'DisplayName', algos{a}.name);
                end
            end
            hold off; grid on;
            
            if m == 1
                switch att_type
                    case 'contiguous_crop', pretty = 'Wycięcie Ciągłe';
                    case 'crop',            pretty = 'Wycięcie Losowe';
                    case 'tamper',          pretty = 'Podmiana (Tamper)';
                    case 'salt_pepper',     pretty = 'Szum Sól/Pieprz';
                    case 'gaussian',        pretty = 'Szum Gaussa';
                    otherwise,              pretty = att_type;
                end
                title(pretty, 'FontWeight', 'bold');
            end
            if t == 1, ylabel(titles{m}, 'FontWeight', 'bold'); end
            if m == rows, xlabel('Intensywność'); end
            if t == cols && m == 1, legend('Location', 'bestoutside', 'Interpreter', 'none'); end
        end
    end
    
    fname = fullfile(save_dir, ['alteration_metrics_plot_' timestamp '.png']);
    saveas(f, fname);
    fprintf('  -> Wykres metryk zapisany: %s\n', fname);
end

function plot_attack_examples(visual_data, algos, attack_configs, save_dir, timestamp)
    num_attacks = size(attack_configs, 1);
    
    for a = 1:length(algos)
        f = figure('Name', ['Wizualizacja: ' algos{a}.name], 'Color', 'w', 'Visible', 'off', 'Position', [100 100 1200 400]);
        
        has_data = false;
        for t = 1:num_attacks, if ~isempty(visual_data{a, t}), has_data = true; break; end; end
        
        if has_data
            for t = 1:num_attacks
                data_struct = visual_data{a, t};
                if isempty(data_struct), continue; end
                
                att_type = attack_configs{t, 1};
                switch att_type
                    case 'contiguous_crop', att_name = 'Ciągłe';
                    case 'crop',            att_name = 'Losowe';
                    case 'tamper',          att_name = 'Podmiana';
                    case 'salt_pepper',     att_name = 'Sól/Pieprz';
                    case 'gaussian',        att_name = 'Gauss';
                    otherwise,              att_name = att_type;
                end
                
                % Ciphertext
                subplot(2, num_attacks, t);
                imshow(data_struct.C_img);
                title(sprintf('%s (%.2f)\nSzyfrogram', att_name, data_struct.intensity), 'FontSize', 8);
                
                % Decrypted
                subplot(2, num_attacks, t + num_attacks);
                imshow(data_struct.P_img);
                title('Odtworzony', 'FontSize', 8);
            end
            sgtitle(['Algorytm: ' algos{a}.name], 'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
            
            safe_name = regexprep(algos{a}.name, '[^a-zA-Z0-9]', '_');
            fname = fullfile(save_dir, ['vis_' safe_name '_' timestamp '.png']);
            saveas(f, fname);
            close(f);
        else
            close(f);
        end
    end
    fprintf('  -> Zapisano wykresy wizualizacji dla każdego algorytmu.\n');
end