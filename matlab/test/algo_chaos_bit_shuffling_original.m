function [t_enc, t_dec, C1, C2, PT] = algo_chaos_bit_shuffling_original(data)
    % ALGO_CHAOS_BIT_SHUFFLING_ORIGINAL - Wrapper for the original version
    
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    
    % 1. Encryption (Original)
    tic;
    [ct1, keys] = chaos_bit_shuffling('encrypt', img_orig);
    t_enc = toc;
    
    % 2. Encryption for metrics
    [ct2, ~] = chaos_bit_shuffling('encrypt', img_mod);
    
    % 3. Decryption
    tic;
    pt_img = chaos_bit_shuffling('decrypt', ct1, keys);
    t_dec = toc;
    
    % 4. Result conversion
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_img);
end
