function [t_enc, t_dec, C1, C2, PT] = algo_chaos_matlab(data)
    % Circular Chaotic Map (Pure MATLAB)
    
    % --- OPTIMIZATION: Extract data BEFORE tic ---
    % Accessing structure fields like data.img_orig inside tic can skew timing
    input_orig = data.img_orig;
    input_mod  = data.img_mod;
    
    tic;
    % Use local variable input_orig
    [C1, k1] = encrypt_circular_chaotic(input_orig);
    t_enc = toc;
    
    % Second run (metrics only, using local variable)
    [C2, ~]  = encrypt_circular_chaotic(input_mod);
    
    tic;
    PT = decrypt_circular_chaotic(C1, k1);
    t_dec = toc;
end
