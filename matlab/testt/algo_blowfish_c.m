function [t_enc, t_dec, C1, C2, PT] = algo_blowfish_c(data)
    % Blowfish-CBC (Software/C implementation via MEX)
    
    key = uint8(randi([0 255], 1, 16));
    iv  = uint8(randi([0 255], 1, 8)); 
    
    tic;
    ct1_flat = crypto_mex(data.img_flat, key, iv, 'bf-cbc', true);
    t_enc = toc;
    
    ct2_flat = crypto_mex(data.img_flat_mod, key, iv, 'bf-cbc', true);
    
    tic;
    pt_flat = crypto_mex(ct1_flat, key, iv, 'bf-cbc', false);
    t_dec = toc;
    
    sz = data.size;
    C1 = reshape(ct1_flat(1:prod(sz)), sz);
    C2 = reshape(ct2_flat(1:prod(sz)), sz);
    PT = reshape(pt_flat(1:prod(sz)), sz);
end
