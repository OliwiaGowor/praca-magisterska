function [t_enc, t_dec, C1, C2, PT] = algo_aes_pixel(data)
    % AES in "Pixel" mode (ECB)
    
    % Generate 256-bit key
    key = uint8(randi([0 255], 1, 32)); 
    
    input_flat = data.img_flat;     % 1D vector
    input_mod  = data.img_flat_mod; % 1D vector
    sz = data.size;
    
    % Array initialization (done once, outside encryption timing per se,
    % or included - depends on methodology. Here included for fairness against DES).
    [w, Nr, s_box, inv_s_box, m9, m11, m13, m14] = aes_init(key);
    
    % --- ENCRYPTION (Vectorized) ---
    tic;
    % Padding the whole vector once.
    [padded_input, pad_len] = pkcs7_pad(input_flat);
    
    ct1_padded = aes_cipher_vectorized(padded_input, w, Nr, s_box);
    t_enc = toc;
    
    % Encrypting the second image
    [padded_mod, ~] = pkcs7_pad(input_mod);
    ct2_padded = aes_cipher_vectorized(padded_mod, w, Nr, s_box);
    
    % --- DECRYPTION ---
    tic;
    pt_padded = aes_inv_cipher_vectorized(ct1_padded, w, Nr, inv_s_box, m9, m11, m13, m14);
    
    % Removing padding
    if pad_len > 0
        pt_flat = pt_padded(1:end-pad_len);
    else
        pt_flat = pt_padded;
    end
    t_dec = toc;
    
    % --- Result formatting ---
    % Cropping ciphertext to image size for visualization purposes
    num_pixels = prod(sz);
    C1 = reshape(ct1_padded(1:num_pixels), sz);
    C2 = reshape(ct2_padded(1:num_pixels), sz);
    
    % Decrypted image
    PT = reshape(pt_flat(1:num_pixels), sz);
end

function [padded, p_len] = pkcs7_pad(data)
    % Helper function for padding
    len = length(data);
    p_len = 16 - mod(len, 16);
    if p_len == 0, p_len = 16; end
    padded = [data, repmat(uint8(p_len), 1, p_len)];
end
