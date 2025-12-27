function [t_enc, t_dec, C1, C2, PT] = algo_hyperchaotic_2d_matlab(data)
    % ALGO_HYPERCHAOTIC_2D - Benchmark wrapper for 2D Chaos Algorithm
    % Compliant with test standard (algo_*.m)
    
    % 1. Prepare Data
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    
    % Algorithm works on double/grayscale, adapter handles it.
    
    % 2. Encryption (Original) - Measure time and retrieve keys
    tic;
    [ct1, keys] = chaos_2d_adapter('encrypt', img_orig);
    t_enc = toc;
    
    % 3. Encryption (Modified) - Only for NPCR/UACI metrics
    % Keys are ignored as this is only for ciphertext comparison
    [ct2, ~] = chaos_2d_adapter('encrypt', img_mod);
    
    % 4. Decryption
    tic;
    pt_img = chaos_2d_adapter('decrypt', ct1, keys);
    t_dec = toc;
    
    % 5. Format Results (uint8 required for metrics)
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_img);
end
