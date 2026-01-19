function output = aes_mode(data, key, iv, mode, algorithm)
% AES_MODE Encrypt/Decrypt using CBC or CTR mode
%
%   output = aes_mode(data, key, iv, mode, algorithm)
%
%   Input:
%       data:      uint8 vector or string
%       key:       AES Key (16, 24, or 32 bytes)
%       iv:        IV (16 bytes)
%       mode:      'encrypt' or 'decrypt'
%       algorithm: 'cbc' (default) or 'ctr'

    % --- Init ---
    if nargin < 5, algorithm = 'cbc'; end
    
    if ischar(data) || isstring(data), data = uint8(char(data)); end
    key = uint8(key);
    iv  = uint8(iv);

    % Load Schedule & Tables once
    [w, Nr, s_box, inv_s_box, m9, m11, m13, m14] = aes_init(key); 

    % --- Processing ---
    if strcmp(algorithm, 'cbc')
        if strcmp(mode, 'encrypt')
            output = cbc_encrypt_loop(data, w, Nr, s_box, iv);
        elseif strcmp(mode, 'decrypt')
            % Use the optimized vectorized decryption
            output = cbc_decrypt_vectorized(data, w, Nr, inv_s_box, iv, m9, m11, m13, m14);
        end
        
    elseif strcmp(algorithm, 'ctr')
        % CTR mode handles encryption and decryption identically
        % It turns AES into a stream cipher and DOES NOT require padding.
        output = ctr_process_loop(data, w, Nr, s_box, iv);
        
    else
        error('Nieznany algorytm. Użyj "cbc" lub "ctr".');
    end
end

%% --- CBC Helper Functions ---

function ciphertext = cbc_encrypt_loop(plaintext, w, Nr, s_box, iv)
    % 1. PKCS#7 Padding
    pad_len = 16 - mod(length(plaintext), 16);
    if pad_len == 0, pad_len = 16; end
    padded_data = [plaintext, repmat(uint8(pad_len), 1, pad_len)];
    
    % 2. Sequential Loop (Block N depends on N-1)
    num_blocks = length(padded_data) / 16;
    ciphertext = zeros(1, length(padded_data), 'uint8');
    previous_block = iv;
    
    for i = 1:num_blocks
        idx_s = (i-1)*16 + 1;
        idx_e = i*16;
        
        input_block = bitxor(padded_data(idx_s:idx_e), previous_block);
        encrypted_block = aes_cipher_core(input_block, w, Nr, s_box);
        
        ciphertext(idx_s:idx_e) = encrypted_block;
        previous_block = encrypted_block;
    end
end

function plaintext = cbc_decrypt_vectorized(ciphertext, w, Nr, inv_s_box, iv, m9, m11, m13, m14)
    if mod(length(ciphertext), 16) ~= 0
        error('Długość szyfrogramu musi być wielokrotnością 16 dla deszyfrowania CBC.');
    end

    % 1. Batch Decrypt
    % Decrypts all blocks simultaneously using the vectorized core
    decrypted_stream = aes_inv_cipher_vectorized(ciphertext, w, Nr, inv_s_box, m9, m11, m13, m14);
    
    % 2. CBC XOR Step (Vectorized)
    % The "previous block" for block N is ciphertext block N-1
    prev_blocks = [iv, ciphertext(1:end-16)];
    padded_plaintext = bitxor(decrypted_stream, prev_blocks);
    
    % 3. Validate & Remove PKCS#7 Padding
    pad_len = double(padded_plaintext(end));
    
    if pad_len > 16 || pad_len == 0
        warning('Nieprawidłowa długość paddingu. Wynik prawdopodobnie uszkodzony.');
        plaintext = padded_plaintext; return;
    end
    
    padding_bytes = padded_plaintext(end-pad_len+1 : end);
    if any(padding_bytes ~= pad_len)
        warning('Nieprawidłowe bajty paddingu. Wynik prawdopodobnie uszkodzony.');
        plaintext = padded_plaintext; return;
    end
    
    plaintext = padded_plaintext(1 : end-pad_len);
end

%% --- CTR Helper Functions ---

function output = ctr_process_loop(data, w, Nr, s_box, iv)
    % CTR Mode: Stream Cipher behavior. 
    % Note: CTR always uses the ENCRYPTION core (aes_cipher_core), never inverse.
    
    num_blocks = ceil(length(data) / 16);
    output = zeros(1, length(data), 'uint8');
    counter_block = iv; 
    
    for i = 1:num_blocks
        idx_start = (i-1)*16 + 1;
        idx_end   = min(i*16, length(data));
        
        % 1. Generate Keystream (Encrypt the Counter)
        keystream = aes_cipher_core(counter_block, w, Nr, s_box);
        
        % 2. XOR keystream with data
        len = idx_end - idx_start + 1;
        input_chunk = data(idx_start:idx_end);
        output_chunk = bitxor(input_chunk, keystream(1:len));
        
        output(idx_start:idx_end) = output_chunk;
        
        % 3. Increment Counter (Big Endian 128-bit integer)
        for b = 16:-1:1
            val = double(counter_block(b)) + 1;
            if val <= 255
                counter_block(b) = uint8(val);
                break; 
            else
                counter_block(b) = 0; % Carry over
            end
        end
    end
end