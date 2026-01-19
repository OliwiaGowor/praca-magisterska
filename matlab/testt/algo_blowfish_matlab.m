function [t_enc, t_dec, C1, C2, PT] = algo_blowfish_matlab(data)
    % ALGO_BLOWFISH_MATLAB Benchmark for Blowfish (Pure MATLAB)
    
    % 1. Setup
    % Blowfish supports variable key lengths (e.g., 16 bytes = 128 bits)
    key = uint8(randi([0 255], 1, 16)); 
    % Blowfish has a 64-bit block size, so IV must be 8 bytes
    iv  = uint8(randi([0 255], 1, 8));  
    
    % Initialize engine (using the ImageBlowfish wrapper)
    bf = ImageBlowfish(key);
    
    % Prepare data (ImageBlowfish accepts image matrix)
    img_orig = data.img_orig; 
    img_mod  = data.img_mod;
    sz       = data.size;
    num_pixels = prod(sz);

    % 2. ENCRYPTION Benchmark
    tic;
    % encryptImage returns a flat byte vector (with padding)
    [ct1_flat, orig_sz, orig_cls] = bf.encryptImage(img_orig, iv);
    t_enc = toc;
    
    % Encrypt modified image (without timing, for correlation/NPCR analysis)
    [ct2_flat, ~, ~] = bf.encryptImage(img_mod, iv);
    
    % 3. DECRYPTION Benchmark
    tic;
    % ImageBlowfish.decryptImage handles reshaping to original size
    PT = bf.decryptImage(ct1_flat, iv, orig_sz, orig_cls);
    t_dec = toc;
    
    % 4. Format outputs for visualization
    % Ciphertext (ct1_flat) has padding, so it is longer than the image.
    % Crop it to the image size only for imshow() visualization purposes.
    
    limit_c1 = min(length(ct1_flat), num_pixels);
    C1 = reshape(ct1_flat(1:limit_c1), sz);
    
    limit_c2 = min(length(ct2_flat), num_pixels);
    C2 = reshape(ct2_flat(1:limit_c2), sz);
    
    % PT is already a matrix thanks to the decryptImage method
end
