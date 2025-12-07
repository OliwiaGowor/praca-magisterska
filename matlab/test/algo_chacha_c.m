function [t_enc, t_dec, C1, C2, PT] = algo_chacha_c(data)
    % ChaCha20 (Software/C implementation via MEX)
    
    key = uint8(randi([0 255], 1, 32));
    iv  = uint8(randi([0 255], 1, 16));
    
    tic;
    ct1_flat = crypto_mex(data.img_flat, key, iv, 'chacha20', true);
    t_enc = toc;
    
    ct2_flat = crypto_mex(data.img_flat_mod, key, iv, 'chacha20', true);
    
    tic;
    % ChaCha20 is symmetric; encryption function also decrypts
    pt_flat = crypto_mex(ct1_flat, key, iv, 'chacha20', true);
    t_dec = toc;
    
    sz = data.size;
    C1 = reshape(ct1_flat, sz);
    C2 = reshape(ct2_flat, sz);
    PT = reshape(pt_flat, sz);
end
