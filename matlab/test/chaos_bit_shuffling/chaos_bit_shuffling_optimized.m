function [out1, out2] = chaos_bit_shuffilng_optimized(mode, varargin)
    % DYSPOZYTOR: Pozwala wywołać funkcje lokalne z zewnątrz
    if strcmp(mode, 'encrypt')
        [out1, out2] = encrypt_image_opt(varargin{1});
    elseif strcmp(mode, 'decrypt')
        out1 = decrypt_image_opt(varargin{1}, varargin{2});
        out2 = []; 
    else
        error('Tryb musi byc "encrypt" lub "decrypt"');
    end
end

% =========================================================
% ENCRYPTION
% =========================================================
function [Encryptedim, keys] = encrypt_image_opt(input_image)
    % --- Preprocessing ---
    if size(input_image, 3) == 3
        A = double(rgb2gray(input_image));
    else
        A = double(input_image);
    end
    
    [rows, cols] = size(A);
    
    % Bit decomposition
    Abits = zeros(rows, cols, 8);
    for b = 1:8
        Abits(:,:,b) = bitget(A, b);
    end
    
    shuffled = Abits; 
    
    % --- Key Generation / Initialization ---
    total_sum = sum(A(:));
    AVG = mod(total_sum/(rows*cols), 1);
    
    if AVG == 0
        AVG = mod(log(rows*cols + total_sum/(rows*cols)), 1);
    end
    
    x_val = AVG;                mu1 = 900 + AVG;    a1 = 2*pi*AVG;
    y_val = mod(10^6*AVG, 1);   mu2 = 901 + y_val;  a2 = 2*pi*y_val;
    z_val = mod(10^9*AVG, 1);   mu3 = 902 + z_val;  a3 = 2*pi*z_val;
    v_val = cos(AVG);           mu4 = 905 + v_val;  a4 = 2*pi*v_val;

    keys = [x_val, mu1, a1, y_val, mu2, a2, z_val, mu3, a3, v_val, mu4, a4];
    
    % --- Shuffle Rows ---
    p_cols = zeros(1, cols);
    p_bits = zeros(1, 8);
    curr_x = x_val; 
    
    for R = 1:rows
        if R ~= 1
            curr_x = cos(mu1*(curr_x^3+curr_x)+a1);
        end
        [p_cols, curr_x] = generate_permutation(curr_x, mu1, a1, cols);
        [p_bits, curr_x] = generate_permutation(curr_x, mu1, a1, 8);
        
        slice = shuffled(R, :, :); 
        slice = slice(1, p_cols, p_bits);
        shuffled(R, :, :) = slice;
    end
    
    % --- Shuffle Columns ---
    curr_y = y_val;
    for C = 1:cols
        if C ~= 1
            curr_y = cos(mu2*(curr_y^3+curr_y)+a2);
        end
        [p_rows, curr_y] = generate_permutation(curr_y, mu2, a2, rows);
        [p_bits, curr_y] = generate_permutation(curr_y, mu2, a2, 8);
        
        slice = shuffled(:, C, :);
        slice = slice(p_rows, 1, p_bits);
        shuffled(:, C, :) = slice;
    end
    
    % --- Shuffle Pixel Level (Depth) ---
    curr_z = z_val;
    for P = 1:8
        if P ~= 1
            curr_z = cos(mu3*(curr_z^3+curr_z)+a3);
        end
        [p_rows, curr_z] = generate_permutation(curr_z, mu3, a3, rows);
        [p_cols, curr_z] = generate_permutation(curr_z, mu3, a3, cols);
        
        slice = shuffled(:, :, P);
        slice = slice(p_rows, p_cols);
        shuffled(:, :, P) = slice;
    end
    
    % --- Diffusion (XOR) ---
    shuffledstream = reshape(shuffled, 1, []);
    num_elements = rows*cols*8;
    stream = false(1, num_elements);
    stream(1) = 1; 
    
    temp_v = v_val;
    mu4_v = mu4; a4_v = a4;
    
    for i = 2:num_elements
        next_v = cos(mu4_v * (temp_v^3 + temp_v) + a4_v);
        val = 10^12 * abs(next_v + temp_v);
        val = val - 2*floor(val/2); 
        stream(i) = floor(val);
        temp_v = next_v;
    end
    
    encryptedstream = xor(shuffledstream, stream);
    Encryptedim = rebuild_image_from_bits(encryptedstream, rows, cols);
end

% =========================================================
% DECRYPTION
% =========================================================
function Reconstructedim = decrypt_image_opt(Encryptedim, keys)
    Encryptedim = double(Encryptedim);
    [rows, cols] = size(Encryptedim);
    
    Encryptedimbits = zeros(rows, cols, 8);
    for b = 1:8
        Encryptedimbits(:,:,b) = bitget(Encryptedim, b);
    end
    encryptedstream = reshape(Encryptedimbits, 1, []);
    
    x_val = keys(1); mu1 = keys(2); a1 = keys(3);
    y_val = keys(4); mu2 = keys(5); a2 = keys(6);
    z_val = keys(7); mu3 = keys(8); a3 = keys(9);
    v_val = keys(10); mu4 = keys(11); a4 = keys(12);
    
    % --- Reverse Diffusion ---
    num_elements = rows*cols*8;
    stream = false(1, num_elements);
    stream(1) = 1;
    temp_v = v_val;
    
    for i = 2:num_elements
        next_v = cos(mu4 * (temp_v^3 + temp_v) + a4);
        val = 10^12 * abs(next_v + temp_v);
        val = val - 2*floor(val/2);
        stream(i) = floor(val);
        temp_v = next_v;
    end
    
    Dencryptedstream = xor(encryptedstream, stream);
    Dshuffledim = reshape(Dencryptedstream, rows, cols, 8);
    
    % --- Reverse Shuffle: Pixel Level ---
    curr_z = z_val;
    for P = 1:8
        if P ~= 1
            curr_z = cos(mu3*(curr_z^3+curr_z)+a3);
        end
        [p_rows, curr_z] = generate_permutation(curr_z, mu3, a3, rows);
        [p_cols, curr_z] = generate_permutation(curr_z, mu3, a3, cols);
        
        slice = Dshuffledim(:, :, P);
        temp_slice = zeros(rows, cols);
        temp_slice(p_rows, p_cols) = slice;
        Dshuffledim(:, :, P) = temp_slice;
    end
    
    % --- Reverse Shuffle: Columns ---
    curr_y = y_val;
    for C = 1:cols
        if C ~= 1
            curr_y = cos(mu2*(curr_y^3+curr_y)+a2);
        end
        [p_rows, curr_y] = generate_permutation(curr_y, mu2, a2, rows);
        [p_bits, curr_y] = generate_permutation(curr_y, mu2, a2, 8);
        
        slice = Dshuffledim(:, C, :);
        temp_slice = zeros(rows, 1, 8);
        temp_slice(p_rows, 1, p_bits) = slice;
        Dshuffledim(:, C, :) = temp_slice;
    end
    
    % --- Reverse Shuffle: Rows ---
    curr_x = x_val;
    for R = 1:rows
        if R ~= 1
            curr_x = cos(mu1*(curr_x^3+curr_x)+a1);
        end
        [p_cols, curr_x] = generate_permutation(curr_x, mu1, a1, cols);
        [p_bits, curr_x] = generate_permutation(curr_x, mu1, a1, 8);
        
        slice = Dshuffledim(R, :, :);
        temp_slice = zeros(1, cols, 8);
        temp_slice(1, p_cols, p_bits) = slice;
        Dshuffledim(R, :, :) = temp_slice;
    end
    
    Reconstructedim = rebuild_image_from_bits(reshape(Dshuffledim, 1, []), rows, cols);
    Reconstructedim = uint8(Reconstructedim);
end

% =========================================================
% FUNKCJE POMOCNICZE
% =========================================================

function [perm, last_val] = generate_permutation(start_val, mu, a, n)
    perm = zeros(1, n);
    used = false(1, n); 
    curr = start_val;
    
    % --- FIX: Zabezpieczenie przed wyjściem poza zakres [1..n] ---
    pos = floor(n * abs(curr)) + 1;
    if pos > n, pos = n; end % Safety clamp
    
    perm(1) = pos;
    used(pos) = true;
    count = 1;
    
    while count < n
        curr = cos(mu * (curr^3 + curr) + a);
        pos = floor(n * abs(curr)) + 1;
        
        % --- FIX: Safety clamp ---
        if pos > n, pos = n; end
        
        if ~used(pos)
            count = count + 1;
            perm(count) = pos;
            used(pos) = true;
        end
    end
    last_val = curr;
end

function img = rebuild_image_from_bits(stream, rows, cols)
    bits_matrix = reshape(stream, rows, cols, 8);
    weights = 2 .^ (0:7); 
    weights = reshape(weights, 1, 1, 8);
    img = sum(bits_matrix .* weights, 3);
end