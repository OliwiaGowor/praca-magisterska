function [t_enc, t_dec, C1, C2, PT] = algo_hu_tian(data)
    % ALGO_HU_TIAN - Wrapper for Hu & Tian Algorithm (Two-Stage Logistic + M-Sequence)
    
    % 1. Prepare Data
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    
    % 2. Encryption
    tic;
    [ct1, keys] = hu_tian_adapter('encrypt', img_orig);
    t_enc = toc;
    
    % 3. Encryption for metrics
    [ct2, ~] = hu_tian_adapter('encrypt', img_mod);
    
    % 4. Decryption - Timing
    tic;
    pt_img = hu_tian_adapter('decrypt', ct1, keys);
    t_dec = toc;
    
    % 5. Format Results
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_img);
end
