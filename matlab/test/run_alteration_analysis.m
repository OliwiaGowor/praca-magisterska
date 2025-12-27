function run_alteration_analysis(data)
    % RUN_ALTERATION_ANALYSIS Comprehensive Alteration Attack Analysis (All Algorithms).
    %
    % Tests robustness of 12 algorithms against:
    % - Contiguous Crop (Simulation of data loss/packet loss)
    % - Random Crop (Noise-like data loss)
    % - Tampering (Substitution)
    % - Noise (Salt & Pepper, Gaussian)
    %
    % Generates:
    % 1. Metrics plots (PSNR, SSIM, SD) vs Intensity.
    % 2. Visual examples (Corrupted Ciphertext vs Reconstructed Image).
    
        fprintf('\n=== ROZPOCZYNAM KOMPLEKSOWĄ ANALIZĘ ATAKÓW (12 ALGORYTMÓW) ===\n');
        
        % --- Suppress known warnings during attacks ---
        orig_warn_state_aes = warning('off', 'AES:BadPadding');
        orig_warn_state_des = warning('off', 'MATLAB:subsassigndimmismatch'); % Occasional shape mismatch in DES
        cleanupObj = onCleanup(@() restore_warnings(orig_warn_state_aes, orig_warn_state_des)); 
    
        % --- 1. Global Configuration ---
        target_size = data.size;
        target_numel = prod(target_size);
        
        % Generate keys for standard algorithms
        key_32 = uint8(randi([0 255], 1, 32)); 
        key_16 = uint8(randi([0 255], 1, 16));
        iv_16  = uint8(randi([0 255], 1, 16));
        iv_8   = uint8(randi([0 255], 1, 8));
        iv_12  = uint8(randi([0 255], 1, 12)); % Nonce for ChaCha
        
        % DES Keys (Hex strings)
        key_des = '1234567890ABCDEF'; 
        iv_des  = 'FEDCBA0987654321';
    
        % --- 2. Algorithm Definitions ---
        algos = {};
        
        % --- GROUP 1: C/MEX Implementations ---
        algos{end+1} = struct('name', 'AES-CBC (MEX)', ...
            'enc', @(img) crypto_mex(img(:)', key_32, iv_16, 'aes-256-cbc', true), ...
            'dec', @(ct)  decrypt_and_restore(ct, @(c) crypto_mex(c, key_32, iv_16, 'aes-256-cbc', false), target_size));
    
        algos{end+1} = struct('name', 'Blowfish (MEX)', ...
            'enc', @(img) crypto_mex(img(:)', key_16, iv_8, 'bf-cbc', true), ...
            'dec', @(ct)  decrypt_and_restore(ct, @(c) crypto_mex(c, key_16, iv_8, 'bf-cbc', false), target_size));
    
        algos{end+1} = struct('name', 'ChaCha20 (MEX)', ...
            'enc', @(img) crypto_mex(img(:)', key_32, iv_12, 'chacha20', true), ...
            'dec', @(ct)  decrypt_and_restore(ct, @(c) crypto_mex(c, key_32, iv_12, 'chacha20', false), target_size));
    
        % --- GROUP 2: MATLAB Standard Implementations ---
        algos{end+1} = struct('name', 'AES-CBC (MATLAB)', ...
            'enc', @(img) aes_mode(img(:)', key_32, iv_16, 'encrypt', 'cbc'), ...
            'dec', @(ct)  decrypt_and_restore(ct, @(c) aes_mode(c, key_32, iv_16, 'decrypt', 'cbc'), target_size));
    
        % Blowfish (Matlab) - Uses ImageBlowfish class
        bf_engine = ImageBlowfish(key_16);
        algos{end+1} = struct('name', 'Blowfish (MATLAB)', ...
            'enc', @(img) bf_engine.encryptImage(img, iv_8), ... % Returns flat vector
            'dec', @(ct)  bf_engine.decryptImage(ct, iv_8, target_size, 'uint8'));
    
        algos{end+1} = struct('name', 'ChaCha20 (MATLAB)', ...
            'enc', @(img) chacha20(key_32, iv_12, 0, img(:)), ...
            'dec', @(ct)  reshape(chacha20(key_32, iv_12, 0, ct), target_size));
    
        % DES (Matlab) - Optimized class
        % Note: DES encrypt returns [ct_bytes, C_img]. We use ct_bytes (vector) for attacks.
        algos{end+1} = struct('name', 'DES-CBC (MATLAB)', ...
            'enc', @(img) des_enc_helper(img, key_des, iv_des), ...
            'dec', @(ct)  DES_ImgProcess_Opt.decrypt(ct, key_des, iv_des, target_size));
    
        % --- GROUP 3: Chaotic Implementations ---
        
        % Chaos Circular (Matlab)
        algos{end+1} = struct('name', 'Chaos Circular', 'type', 'chaos_stateful', ...
            'enc_gen', @(img) encrypt_circular_chaotic(img), ... % Returns [C, keys]
            'dec_gen', @(ct, keys) decrypt_circular_chaotic(ct, keys));
    
        % Chaos Bit Shuffling (Moysis) - FIXED: Wraps output in uint8()
        algos{end+1} = struct('name', 'Chaos Bit Shuffling', 'type', 'chaos_stateful', ...
            'enc_gen', @(img) bit_shuffling_enc_wrapper(img), ... 
            'dec_gen', @(ct, keys) chaos_bit_shuffling_optimized('decrypt', ct, keys));
            
        % Hyperchaotic 2D (Liu)
        algos{end+1} = struct('name', 'Hyperchaotic 2D', 'type', 'chaos_stateful', ...
            'enc_gen', @(img) chaos_2d_adapter('encrypt', img), ...
            'dec_gen', @(ct, keys) chaos_2d_adapter('decrypt', ct, keys));
    
        % Hu & Tian (Two-Stage Logistic)
        algos{end+1} = struct('name', 'Hu & Tian', 'type', 'chaos_stateful', ...
            'enc_gen', @(img) hu_tian_adapter('encrypt', img), ...
            'dec_gen', @(ct, keys) hu_tian_adapter('decrypt', ct, keys));
    
        % Chaos Maps (3D Logistic + Chirikov)
        % Requires local helper functions (defined at bottom)
        algos{end+1} = struct('name', 'Chaos Maps (3D Log)', 'type', 'chaos_stateful', ...
            'enc_gen', @(img) chaos_map_enc_helper(img), ...
            'dec_gen', @(ct, state) chaos_map_dec_helper(ct, state, target_size));
    
    
        % --- 3. Attack Configuration ---
        attack_configs = {
            'contiguous_crop', [0.1, 0.2, 0.3, 0.5];    % "Black Hole" attack
            'crop',            [0.1, 0.2, 0.3, 0.5];    % Noise-like crop
            'tamper',          [0.05, 0.1, 0.2, 0.4];   % Image substitution
            'salt_pepper',     [0.01, 0.05, 0.1, 0.2];  % Salt & Pepper noise
            'gaussian',        [0.005, 0.01, 0.02, 0.05]; % Gaussian noise
        };
        
        vis_intensities = containers.Map;
        vis_intensities('contiguous_crop') = 0.2;
        vis_intensities('crop')            = 0.2;
        vis_intensities('tamper')          = 0.1;
        vis_intensities('salt_pepper')     = 0.05;
        vis_intensities('gaussian')        = 0.01;
    
        % Load intruder image
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
            fprintf('[%d/%d] Algorytm: %-25s ... ', a, length(algos), algo.name);
            
            % 4a. Initial Encryption (Get Ciphertext and State/Keys)
            try
                if isfield(algo, 'type') && strcmp(algo.type, 'chaos_stateful')
                    % Stateful algorithms return [Ciphertext, State]
                    [C_base, state] = algo.enc_gen(data.img_orig);
                    decrypt_fn = @(ct) algo.dec_gen(ct, state);
                else
                    % Stateless/Standard algorithms
                    C_base = algo.enc(data.img_orig);
                    decrypt_fn = algo.dec;
                end
            catch ME
                fprintf('NIEPOWODZENIE (Szyfr): %s\n', ME.message);
                continue;
            end
            
            % 4b. Attack Loop
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
                    
                    % Apply Attack
                    C_prim = apply_attack(C_base, att_type, inte, img_tamper);
                    
                    % Decrypt
                    try
                        P_prim = decrypt_fn(C_prim);
                    catch
                        % Decryption failed (padding error etc.) -> return black image
                        P_prim = zeros(target_size, 'uint8');
                    end
                    
                    % Shape correction
                    if ~isequal(size(P_prim), target_size)
                        P_tmp = zeros(target_size, 'uint8');
                        num_copy = min(numel(P_prim), target_numel);
                        P_prim_flat = P_prim(:);
                        P_tmp(1:num_copy) = P_prim_flat(1:num_copy);
                        P_prim = P_tmp;
                    end
                    
                    % Calculate Metrics
                    [p, s, d] = calculate_alteration_metrics(data.img_orig, P_prim);
                    psnr_vec(i) = p; ssim_vec(i) = s; sd_vec(i) = d;
                    
                    % Store for visualization
                    if abs(inte - vis_int) < 0.001 || i == 1
                         % Convert C_prim to 2D for display if it is a vector
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
                
                % Store Results
                key_res = sprintf('a%d_%s', a, att_type);
                results.(key_res).psnr = psnr_vec;
                results.(key_res).ssim = ssim_vec;
                results.(key_res).sd   = sd_vec;
                results.(key_res).intensities = intensities;
            end
            fprintf('Gotowe.\n');
        end
        
        fprintf('Generowanie wykresów...\n');
        plot_alteration_results(results, algos, attack_configs);
        plot_attack_examples(visual_data, algos, attack_configs);
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
        out = ct; % Return vector
    end
    
    % --- Wrapper for Chaos Bit Shuffling ---
    function [C, keys] = bit_shuffling_enc_wrapper(img)
        % chaos_bit_shuffling_optimized returns double [0, 255].
        % We must cast to uint8 for correct display in imshow.
        [C_double, keys] = chaos_bit_shuffling_optimized('encrypt', img);
        C = uint8(C_double);
    end
    
    % --- Logic for Chaos Maps (3D Logistic + Chirikov) ---
    % Recreated here to allow state capture (keys/seqs)
    
    function [ct, state] = chaos_map_enc_helper(img)
        [rows, cols, ~] = size(img);
        num_pixels = rows * cols;
        key_chirikov = rand(1, 4) * 10; % Random key for Chirikov
        
        % 1. Diffusion (Logistic)
        [diffused_img, seqs] = logistic_diffusion_core(img, num_pixels);
        
        % 2. Confusion (Chirikov)
        [ct, ~] = encrypt(diffused_img, key_chirikov);
        
        state.key = key_chirikov;
        state.seqs = seqs;
    end
    
    function img_out = chaos_map_dec_helper(ct, state, sz)
        % 1. Inverse Confusion
        dec_chirikov = decrypt(ct, state.key);
        
        % 2. Inverse Diffusion
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
    
    function plot_alteration_results(results, algos, attack_configs)
        figure('Name', 'Analiza Metryk', 'Color', 'w', 'Position', [50 50 1400 900]);
        metrics = {'psnr', 'ssim', 'sd'};
        titles = {'PSNR (Im wyżej tym lepiej)', 'SSIM (Im bliżej 1 tym lepiej)', 'SD (Im niżej tym lepiej)'};
        rows = 3; cols = size(attack_configs, 1);
        
        % Colors for 12 algos
        colors = parula(length(algos));
        
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
                    % Translate attack name to Polish for display
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
                if t == cols && m == 1, legend('Location', 'bestoutside'); end
            end
        end
    end
    
    function plot_attack_examples(visual_data, algos, attack_configs)
        num_attacks = size(attack_configs, 1);
        
        % Group algorithms by type to reduce number of figures
        % Or just one figure per algorithm as before
        for a = 1:length(algos)
            f = figure('Name', ['Wizualizacja: ' algos{a}.name], 'Color', 'w', 'Position', [100 100 1200 400]);
            % Check if any data exists
            has_data = false;
            for t = 1:num_attacks, if ~isempty(visual_data{a, t}), has_data = true; break; end; end
            if ~has_data, close(f); continue; end
    
            for t = 1:num_attacks
                data_struct = visual_data{a, t};
                if isempty(data_struct), continue; end
                
                % Translate attack name
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
            sgtitle(['Algorytm: ' algos{a}.name], 'FontWeight', 'bold', 'FontSize', 12);
        end
    end