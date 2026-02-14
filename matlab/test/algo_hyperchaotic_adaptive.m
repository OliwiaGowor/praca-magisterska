function [t_enc, t_dec, C1, C2, PT] = algo_hyperchaotic_adaptive(data)
    % ALGO_HYPERCHAOTIC_ADAPTIVE Adapter for 2D-RA algorithm + Adaptive Diffusion
    
    % Adding path to algorithm modules
    currentPath = fileparts(mfilename('fullpath'));
    algoPath = fullfile(currentPath, 'hyperchaotic_adaptive');
    addpath(algoPath);
    
    img_orig = data.img_orig;
    img_mod = data.img_mod;
    sz = data.size;
    
    tic;
    % 1. Key and sequence generation (dependent on image hash)
    [seq_enc, params_enc] = hc_key_gen(img_orig);
    
    % 2. Encryption
    C1 = hc_encrypt(img_orig, seq_enc);
    t_enc = toc;
    
    %% --- Encryption for metrics ---
    [seq_mod, ~] = hc_key_gen(img_mod);
    C2 = hc_encrypt(img_mod, seq_mod);
    
    %% --- Decryption ---
    tic;
    PT = hc_decrypt(C1, seq_enc, sz);
    t_dec = toc;
    
end
