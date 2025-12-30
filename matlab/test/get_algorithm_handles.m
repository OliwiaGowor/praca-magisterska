function [encrypt_func, decrypt_func] = get_algorithm_handles(algo_name)
    % Mapuje nazwę algorytmu na funkcje szyfrujące i deszyfrujące
    
    switch algo_name
        case 'AES-CBC (Hardware-C)'
            encrypt_func = @(img) aes_wrapper(img, 'enc', 'aes-256-cbc');
            decrypt_func = @(ct, key) aes_wrapper(ct, 'dec', 'aes-256-cbc', key);
            
        case 'Blowfish-CBC (Software-C)'
            encrypt_func = @(img) aes_wrapper(img, 'enc', 'bf-cbc');
            decrypt_func = @(ct, key) aes_wrapper(ct, 'dec', 'bf-cbc', key);
            
        case 'ChaCha20 (Software-C)'
             encrypt_func = @(img) aes_wrapper(img, 'enc', 'chacha20');
             decrypt_func = @(ct, key) aes_wrapper(ct, 'dec', 'chacha20', key);

        case 'Chaotyczny (MATLAB)'
             encrypt_func = @encrypt_circular_chaotic;
             decrypt_func = @decrypt_circular_chaotic;
             
        case 'Chaos Bit Shuffling (MATLAB)'
            encrypt_func = @(img) chaos_bit_shuffilng_optimized('encrypt', img);
            decrypt_func = @(ct, key) chaos_bit_shuffilng_optimized('decrypt', ct, key);
            
        case 'Hyperchaotic 2D (MATLAB)'
            encrypt_func = @(img) chaos_2d_adapter('encrypt', img);
            decrypt_func = @(ct, key) chaos_2d_adapter('decrypt', ct, key);
            
        case 'Hu & Tian - Two-Stage Logistic (MATLAB)'
            encrypt_func = @(img) hu_tian_adapter('encrypt', img);
            decrypt_func = @(ct, key) hu_tian_adapter('decrypt', ct, key);
            
        case 'Entropy Wavelet-Chaos (Pourasad et al.)'
            encrypt_func = @(img) entropy_adapter('encrypt', img);
            decrypt_func = @(ct, key) entropy_adapter('decrypt', ct, key);
            
        otherwise
            error('Algorytm "%s" nie jest obsługiwany w trybie interaktywnym.', algo_name);
    end
end

% Wrapper dla funkcji MEX (C)
function [res, keys] = aes_wrapper(input, mode, algo_type, keys_in)
    input_flat = input(:)';
    sz = size(input);
    
    if strcmp(mode, 'enc')
        if contains(algo_type, 'aes')
            key = uint8(randi([0 255], 1, 32)); iv = uint8(randi([0 255], 1, 16));
        elseif contains(algo_type, 'chacha')
            key = uint8(randi([0 255], 1, 32)); iv = uint8(randi([0 255], 1, 8));
        else % blowfish
            key = uint8(randi([0 255], 1, 16)); iv = uint8(randi([0 255], 1, 8));
        end
        % Wywołanie MEX
        ct = crypto_mex(input_flat, key, iv, algo_type, true);
        res = reshape(ct(1:numel(input)), sz);
        keys = {key, iv};
    else
        key = keys_in{1}; iv = keys_in{2};
        pt = crypto_mex(input_flat, key, iv, algo_type, false);
        res = reshape(pt(1:numel(input)), sz);
        keys = [];
    end
end