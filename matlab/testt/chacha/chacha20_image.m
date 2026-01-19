function [img_encrypted, img_decrypted] = chacha20_image(image_path, key, nonce, counter)
% CHACHA20_IMAGE Robust Image Encryption (Grayscale & RGB)
%
%   [enc, dec] = chacha20_image('photo.jpg', key, nonce, counter)

    % 1. Read Image
    if ischar(image_path) || isstring(image_path)
        img_orig = imread(image_path);
    else
        img_orig = image_path;
    end

    % 2. Store Metadata
    % We need dimensions (Height, Width, Depth) to reconstruct the image later
    dims = size(img_orig);
    
    % 3. Flatten
    % Convert multidimensional matrix to a single 1D byte stream
    flat_input = uint8(img_orig(:));
    
    % 4. Encrypt
    flat_encrypted = chacha20(key, nonce, counter, flat_input);
    
    % 5. Reconstruct Encrypted Image
    img_encrypted = reshape(flat_encrypted, dims);
    
    % 6. Decrypt (Verification Step)
    flat_decrypted = chacha20(key, nonce, counter, flat_encrypted);
    img_decrypted  = reshape(flat_decrypted, dims);
    
    % --- Visualization ---
    figure('Name', 'Szyfrowanie Obrazu ChaCha20');
    
    subplot(1, 3, 1);
    imshow(img_orig);
    title(sprintf('Oryginał (%s)', strjoin(string(dims), 'x')));
    
    subplot(1, 3, 2);
    imshow(img_encrypted);
    title('Zaszyfrowany (Szum)');
    
    subplot(1, 3, 3);
    imshow(img_decrypted);
    title('Odszyfrowany (Odtworzony)');
end
