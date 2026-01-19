function [t_enc, t_dec, C1, C2, PT] = algo_chaos_bit_shuffling(data)
    % Wrapper for Chaos Bit Shuffling (Moysis et al.)
    
    % 1. Data
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    sz = data.size;
    
    % Algorithm works on Grayscale. If input is RGB, 
    % conversion happens inside algorithm, but result will be 2D.
    
    % 2. Encryption
    tic;
    % Call via dispatcher
    [ct1, keys1] = chaos_bit_shuffling_optimized('encrypt', img_orig);
    t_enc = toc;
    
    [ct2, ~] = chaos_bit_shuffling_optimized('encrypt', img_mod);
    
    % 3. Decryption
    tic;
    pt_gray = chaos_bit_shuffling_optimized('decrypt', ct1, keys1);
    t_dec = toc;
    
    % 4. Format
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_gray);
end
