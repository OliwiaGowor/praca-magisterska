function [t_enc, t_dec, C1, C2, PT] = algo_des_c(data)
    % ALGO_DES_C - Wrapper for DES-CBC (OpenSSL via MEX)
    % Requires compiled crypto_mex with 'legacy' support (OpenSSL 3.0+)
    
    % 1. Data preparation (flattening to 1D)
    if isfield(data, 'img_flat')
        p_flat = data.img_flat;
        p_mod_flat = data.img_flat_mod;
    else
        p_flat = data.img_orig(:)';
        p_mod_flat = data.img_mod(:)';
    end
    
    sz = data.size;
    
    key = uint8(randi([0 255], 1, 8));
    iv  = uint8(randi([0 255], 1, 8));
    
    algo_name = 'des-cbc';
    
    % 2. Encryption 
    tic;
    ct1_flat = crypto_mex(p_flat, key, iv, algo_name, true);
    t_enc = toc;
    
    % 3. Encryption (Modified)
    ct2_flat = crypto_mex(p_mod_flat, key, iv, algo_name, true);
    
    % 4. Decryption
    tic;
    pt_flat = crypto_mex(ct1_flat, key, iv, algo_name, false);
    t_dec = toc;
    
    % 5. Result formatting
    % Cropping to image size for visualization and metrics.
    
    len = prod(sz);
    
    % Dimension error protection
    if numel(ct1_flat) >= len
        C1 = reshape(ct1_flat(1:len), sz);
    else
        C1 = zeros(sz, 'uint8');
    end
    
    if numel(ct2_flat) >= len
        C2 = reshape(ct2_flat(1:len), sz);
    else
        C2 = zeros(sz, 'uint8');
    end
    
    if numel(pt_flat) >= len
        PT = reshape(pt_flat(1:len), sz);
    else
        PT = zeros(sz, 'uint8');
    end
end
