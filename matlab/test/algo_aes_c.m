function [t_enc, t_dec, C1, C2, PT] = algo_aes_c(data)
    % AES-256-CBC
    
    key = uint8(randi([0 255], 1, 32)); 
    iv  = uint8(randi([0 255], 1, 16)); 
    
    % Assign data to local variables BEFORE measurement.
    input_plain = data.img_flat;
    input_mod   = data.img_flat_mod;
    
    tic;
    % Pass local variable input_plain instead of data.img_flat
    ct1_flat = crypto_mex(input_plain, key, iv, 'aes-256-cbc', true);
    t_enc = toc;
    
    ct2_flat = crypto_mex(input_mod, key, iv, 'aes-256-cbc', true);
    
    tic;
    pt_flat = crypto_mex(ct1_flat, key, iv, 'aes-256-cbc', false);
    t_dec = toc;
    
    sz = data.size;
    C1 = reshape(ct1_flat(1:prod(sz)), sz);
    C2 = reshape(ct2_flat(1:prod(sz)), sz);
    PT = reshape(pt_flat(1:prod(sz)), sz);
end
