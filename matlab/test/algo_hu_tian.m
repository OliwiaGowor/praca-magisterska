function [t_enc, t_dec, C1, C2, PT] = algo_hu_tian(data)
    % ALGO_HU_TIAN - Wrapper for Hu & Tian Algorithm (Two-Stage Logistic + M-Sequence)
    % Compliant with the test standard.
    
    % 1. Prepare Data
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    
    % 2. Encryption (Original) - Timing
    tic;
    % Adapter returns encrypted image (ct1) and key structure (keys)
    [ct1, keys] = hu_tian_adapter('encrypt', img_orig);
    t_enc = toc;
    
    % 3. Encryption (Modified) - for NPCR/UACI metrics
    % We use the same parameters generated inside the adapter
    % (for valid NPCR tests, the stream generation should be deterministic)
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
