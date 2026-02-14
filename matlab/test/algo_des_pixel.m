function [t_enc, t_dec, C1, C2, PT] = algo_des_pixel(data)
    % DES in "Pixel" mode (ECB)
    
    key_hex = DES_ImgProcess_Opt.generateRandomHex(8);
    
    input_orig = data.img_orig;
    input_mod  = data.img_mod;
    sz = data.size;
    
    % --- ENCRYPTION ---
    tic;
    [ct1_bytes, C1] = encrypt_des_ecb(input_orig, key_hex);
    t_enc = toc;
    
    [ct2_bytes, C2] = encrypt_des_ecb(input_mod, key_hex);
    
    % --- DECRYPTION ---
    tic;
    [~, PT] = decrypt_des_ecb(ct1_bytes, key_hex, sz);
    t_dec = toc;
end

function [cipherBytes, cipherImgVis] = encrypt_des_ecb(imgUint8, keyHex)
    % ECB encryption logic (Pixel-based)
    
    flatPixels = imgUint8(:)'; 
    padLen = 8 - mod(length(flatPixels), 8);
    paddedPixels = [flatPixels, repmat(uint8(padLen), 1, padLen)];
    
    % uint8 -> uint64 (Big Endian)
    blocksIn = swapbytes(typecast(paddedPixels, 'uint64'));
    key64 = DES_ImgProcess_Opt.hexToUint64(keyHex);
    
    % ECB method call
    blocksOut = DES_ImgProcess_Opt.Engine.encrypt_data_ecb(blocksIn, key64);
    
    cipherBytes = typecast(swapbytes(blocksOut), 'uint8');
    
    % Visualization (cropping to original size)
    target_len = numel(imgUint8);
    if length(cipherBytes) >= target_len
        visBytes = cipherBytes(1:target_len);
    else
        visBytes = cipherBytes;
    end
    cipherImgVis = reshape(visBytes, size(imgUint8));
end

function [plainBytes, plainImgVis] = decrypt_des_ecb(cipherBytes, keyHex, originalSize)
    % ECB decryption logic
    
    % uint8 -> uint64 (Big Endian)
    blocksIn = swapbytes(typecast(cipherBytes, 'uint64'));
    key64 = DES_ImgProcess_Opt.hexToUint64(keyHex);
    
    % ECB method call (Decrypt)
    blocksOut = DES_ImgProcess_Opt.Engine.decrypt_data_ecb(blocksIn, key64);
    
    decryptedBytes = typecast(swapbytes(blocksOut), 'uint8');
    
    % Unpad PKCS#7
    try
        padLen = double(decryptedBytes(end));
        if padLen > 0 && padLen <= 8
             % Padding validity check
             padding = decryptedBytes(end-padLen+1:end);
             if all(padding == padLen)
                 plainBytes = decryptedBytes(1:end-padLen);
             else
                 plainBytes = decryptedBytes; % Padding error, return whole
             end
        else
             plainBytes = decryptedBytes;
        end
    catch
        plainBytes = decryptedBytes;
    end
    
    % Image reconstruction
    try
        plainImgVis = reshape(plainBytes, originalSize);
    catch
        % Fallback in case of size error
        warning('Rozmiar po deszyfrowaniu nie pasuje do oryginału.');
        plainImgVis = zeros(originalSize, 'uint8');
        target_len = prod(originalSize);
        if length(plainBytes) >= target_len
            plainImgVis(:) = plainBytes(1:target_len);
        end
    end
end
