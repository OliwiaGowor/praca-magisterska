function [out1, out2] = chaos_bit_shuffling_optimized(mode, varargin)
    % DISPATCHER: Allows calling local functions from outside
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
    
    % Bit decomposition - VECTORIZED
    % Instead of looping over pixels, operate on matrices
    % Abits dimension: rows x cols x 8
    Abits = zeros(rows, cols, 8);
    for b = 1:8
        Abits(:,:,b) = bitget(A, b);
    end
    
    shuffled = Abits; 
    
    % --- Key Generation / Initialization ---
    total_sum = sum(A(:));
    AVG = mod(total_sum/(rows*cols), 1);
    
    % Parameter initialization (as in original)
    if AVG == 0
        AVG = mod(log(rows*cols + total_sum/(rows*cols)), 1);
    end
    
    x_val = AVG;                mu1 = 900 + AVG;    a1 = 2*pi*AVG;
    y_val = mod(10^6*AVG, 1);   mu2 = 901 + y_val;  a2 = 2*pi*y_val;
    z_val = mod(10^9*AVG, 1);   mu3 = 902 + z_val;  a3 = 2*pi*z_val;
    v_val = cos(AVG);           mu4 = 905 + v_val;  a4 = 2*pi*v_val;

    keys = [x_val, mu1, a1, y_val, mu2, a2, z_val, mu3, a3, v_val, mu4, a4];
    
    % Preallocate buffers for chaotic sequences (for speed)
    % We assume pessimistically that the while loop might run more times,
    % but here we allocate working buffers.
    
    % --- Shuffle Rows ---
    % Instead of matrix multiplication: shuffled(R,:,:) = shuffled(R, p_cols, p_bits)
    
    % Buffer for column permutation (different for each row)
    p_cols = zeros(1, cols);
    p_bits = zeros(1, 8);
    
    curr_x = x_val; % Local state copy
    
    for R = 1:rows
        if R ~= 1
            % Update state for new row
            curr_x = cos(mu1*(curr_x^3+curr_x)+a1);
        end
        
        % Generating permutation for columns (cols)
        p_cols = generate_permutation(curr_x, mu1, a1, cols);
        % Update curr_x to the last value from generator
        % (Generator returns permutation and last state, but here
        % we must maintain logic continuity. In original,
        % variable 'x' is overwritten in while loop.
        % Helper 'generate_permutation' returns last x too)
        [p_cols, curr_x] = generate_permutation(curr_x, mu1, a1, cols);
        
        % Generating permutation for bits (8)
        [p_bits, curr_x] = generate_permutation(curr_x, mu1, a1, 8);
        
        % Applying permutation (indexing instead of matrix mult)
        % left * M * right -> M(left_perm, right_perm)
        % Here left is column perm, right is bit perm
        % Note: original was left(cols x cols) * shuffled(R,:,:) * right(8x8)
        % shuffled(R,:,:) is vector 1 x cols x 8.
        % Perm dimensions: (1, cols, 8).
        % left acts on dim 2 (cols), right acts on dim 3 (8).
        
        % Get slice
        slice = shuffled(R, :, :); 
        % Permute
        slice = slice(1, p_cols, p_bits);
        % Save
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
        
        % shuffled(:,C,:) is rows x 1 x 8
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
        
        % shuffled(:,:,P) is rows x cols
        slice = shuffled(:, :, P);
        slice = slice(p_rows, p_cols);
        shuffled(:, :, P) = slice;
    end
    
    % --- Diffusion (XOR) ---
    % Linearization
    shuffledstream = reshape(shuffled, 1, []);
    num_elements = rows*cols*8;
    
    % Diffusion stream generation (preallocation)
    v_seq = zeros(1, num_elements);
    stream = false(1, num_elements); % Logicals are smaller
    
    curr_v = v_val;
    % First element specific (stream(1)=1 in original)
    stream(1) = 1; 
    
    % This is a dependent loop (v(i) depends on v(i-1)), hard to vectorise fully
    % But we can speed up calc by avoiding array resizing
    
    % Start value for loop
    prev_v = curr_v;
    
    % Optimization of diffusion loop
    % Extract constants before loop
    mu4_v = mu4; a4_v = a4;
    
    % Since v(i) depends on v(i-1), we must iterate.
    % We can use v_seq array.
    v_seq(1) = 0; % placeholder, unused in loop from i=2
    
    % Chaos is iterative.
    % v(i) = cos... v(i-1)
    % v(1) is input keys(10). But in loop i goes 2:end.
    % v(1) in loop corresponds to v(i-1) for i=2.
    
    temp_v = prev_v;
    % Using preallocated 'stream' array
    % MATLAB JIT (Just-In-Time) handles preallocated loops well.
    
    for i = 2:num_elements
        % Calculate v(i) based on v(i-1) (temp_v)
        next_v = cos(mu4_v * (temp_v^3 + temp_v) + a4_v);
        
        % Calculate stream bit
        % Original: floor(mod(10^12 * abs(next_v + temp_v), 2))
        val = 10^12 * abs(next_v + temp_v);
        % Faster modulo 2 for positive numbers:
        val = val - 2*floor(val/2); 
        stream(i) = floor(val);
        
        temp_v = next_v;
    end
    
    % XOR
    encryptedstream = xor(shuffledstream, stream);
    
    % --- Reconstruction ---
    % VECTORIZED bit composition
    Encryptedim = rebuild_image_from_bits(encryptedstream, rows, cols);
end

% =========================================================
% DECRYPTION
% =========================================================
function Reconstructedim = decrypt_image_opt(Encryptedim, keys)
    Encryptedim = double(Encryptedim);
    [rows, cols] = size(Encryptedim);
    
    % Bit decomposition - VECTORIZED
    Encryptedimbits = zeros(rows, cols, 8);
    for b = 1:8
        Encryptedimbits(:,:,b) = bitget(Encryptedim, b);
    end
    encryptedstream = reshape(Encryptedimbits, 1, []);
    
    % Parse Keys
    x_val = keys(1); mu1 = keys(2); a1 = keys(3);
    y_val = keys(4); mu2 = keys(5); a2 = keys(6);
    z_val = keys(7); mu3 = keys(8); a3 = keys(9);
    v_val = keys(10); mu4 = keys(11); a4 = keys(12);
    
    % --- Reverse Diffusion (XOR) ---
    % Stream generation identical to encryption
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
        
        % Reversing permutation: Dshuffledim(p_rows, p_cols) = Current
        % Meaning Current(inv_p_rows, inv_p_cols)
        % In MATLAB:
        % If B = A(p), then A(p) = B restores state (assignment by indexing)
        
        slice = Dshuffledim(:, :, P);
        % Here we do reverse: insert values into correct places
        % Assignment method: Temp(p_rows, p_cols) = slice
        
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
    
    % --- Reconstruction ---
    Reconstructedim = rebuild_image_from_bits(reshape(Dshuffledim, 1, []), rows, cols);
    Reconstructedim = uint8(Reconstructedim);
end

% =========================================================
% HELPER FUNCTIONS
% =========================================================

function [perm, last_val] = generate_permutation(start_val, mu, a, n)
    % Generates permutation of numbers 1:n using chaotic map.
    % Optimized version with logical array for O(1) lookups.
    
    perm = zeros(1, n);
    used = false(1, n); % Logical array "is used?"
    
    curr = start_val;
    
    % First element
    pos = floor(n * abs(curr)) + 1;
    perm(1) = pos;
    used(pos) = true;
    
    count = 1;
    
    % While loop until n unique found
    while count < n
        % Chaos map iteration
        curr = cos(mu * (curr^3 + curr) + a);
        pos = floor(n * abs(curr)) + 1;
        
        % Fast check
        if ~used(pos)
            count = count + 1;
            perm(count) = pos;
            used(pos) = true;
        end
    end
    last_val = curr;
end

function img = rebuild_image_from_bits(stream, rows, cols)
    % Fast bit conversion (vectorized bi2de)
    % stream is vector 1 x (rows*cols*8)
    
    % Reshape to matrix: (rows*cols) x 8
    % Note: bitget(A, 1) is LSB.
    % Original code did reshape test(i,j,:) -> 1x1x8.
    % Assuming bit order 1..8 (LSB..MSB).
    
    bits_matrix = reshape(stream, rows, cols, 8);
    
    % Bit weights: 1, 2, 4, 8, 16, 32, 64, 128
    weights = 2 .^ (0:7); % [1 2 4 8 16 32 64 128]
    weights = reshape(weights, 1, 1, 8);
    
    % Weighted sum along 3rd dimension
    img = sum(bits_matrix .* weights, 3);
end
