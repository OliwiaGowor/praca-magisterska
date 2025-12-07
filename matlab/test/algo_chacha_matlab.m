function [t_enc, t_dec, C1, C2, PT] = algo_chacha_matlab(data)
    % ALGO_CHACHA_MATLAB Benchmark wrapper for ChaCha20 (Pure MATLAB)
    
    % 1. Setup
    % ChaCha20 requires a 32-byte key and 12-byte nonce
    key     = uint8(randi([0 255], 1, 32)); 
    nonce   = uint8(randi([0 255], 1, 12));
    counter = 0; 
    
    % Prepare data: flatten image to stream as ChaCha20 works on vectors
    img_orig_flat = data.img_orig(:); 
    img_mod_flat  = data.img_mod(:);
    sz            = data.size;

    % 2. ENCRYPTION Benchmark
    tic;
    ct1_flat = chacha20(key, nonce, counter, img_orig_flat);
    t_enc = toc;
    
    % Encrypt modified image (for NPCR/UACI metrics)
    ct2_flat = chacha20(key, nonce, counter, img_mod_flat);
    
    % 3. DECRYPTION Benchmark
    tic;
    % ChaCha20 is symmetric; applying the same key stream decrypts the data
    pt_flat = chacha20(key, nonce, counter, ct1_flat);
    t_dec = toc;
    
    % 4. Format outputs
    % Reshape flat vectors back to image matrices
    C1 = reshape(ct1_flat, sz);
    C2 = reshape(ct2_flat, sz);
    PT = reshape(pt_flat, sz);
end
