function output = chacha20(key, nonce, counter, input)
% CHACHA20 High-Performance Encrypt/Decrypt (RFC 7539)
%   output = chacha20(key, nonce, counter, input)

    % --- 1. Input Validation and Formatting ---
    if ischar(input) || isstring(input)
        input = uint8(char(input));
    end
    
    % Force column vectors
    key = uint8(key(:));
    nonce = uint8(nonce(:));
    input = uint8(input(:));
    
    if length(key) ~= 32
        error('Key must be 32 bytes.');
    end
    if length(nonce) ~= 12
        error('Nonce must be 12 bytes.');
    end
    
    num_bytes = length(input);
    if num_bytes == 0
        output = uint8([]);
        return;
    end
    
    % Pre-allocate output
    output = zeros(size(input), 'uint8');
    
    % --- 2. Setup Constants and Inputs ---
    consts = uint32([0x61707865; 0x3320646e; 0x79622d32; 0x6b206574]);
    k = typecast_bytes_to_uint32(key);
    n = typecast_bytes_to_uint32(nonce);
    
    % --- 3. Batch Processing Parameters ---
    % Process in chunks to manage memory usage
    BATCH_SIZE = 4096; 
    BLOCK_SIZE = 64;
    num_total_blocks = ceil(num_bytes / BLOCK_SIZE);
    
    % Indices for vectorized Quarter Rounds
    c_a = [1; 2; 3; 4];    c_b = [5; 6; 7; 8];
    c_c = [9; 10; 11; 12]; c_d = [13; 14; 15; 16];
    d_a = [1; 2; 3; 4];    d_b = [6; 7; 8; 5];
    d_c = [11; 12; 9; 10]; d_d = [16; 13; 14; 15];

    % --- 4. Main Processing Loop ---
    for b_start = 1 : BATCH_SIZE : num_total_blocks
        
        % Determine batch range
        b_end = min(b_start + BATCH_SIZE - 1, num_total_blocks);
        current_batch_len = b_end - b_start + 1;
        
        % --- A. Build State Matrix for this Batch ---
        state = zeros(16, current_batch_len, 'uint32');
        state(1:4, :)   = repmat(consts, 1, current_batch_len);
        state(5:12, :)  = repmat(k, 1, current_batch_len);
        state(14:16, :) = repmat(n, 1, current_batch_len);
        
        % Handle Counter Wrapping manually
        % Use double precision for calculation, apply mod 2^32, then cast back.
        rel_indices = 0 : (current_batch_len - 1);
        abs_indices = (b_start - 1) + rel_indices;
        
        % Counter addition with modulo
        state(13, :) = uint32(mod(double(counter) + double(abs_indices), 4294967296)); 

        % --- B. Generate Keystream (Vectorized ARX) ---
        initial_state = state;
        
        for i = 1:10
            state = quarter_round_vect(state, c_a, c_b, c_c, c_d);
            state = quarter_round_vect(state, d_a, d_b, d_c, d_d);
        end
        
        % Final Addition must also be Modulo 2^32
        state = uint32(mod(double(state) + double(initial_state), 4294967296));
        
        % --- C. Serialize and XOR ---
        [~, ~, endian] = computer;
        if endian == 'B'
            state = swapbytes(state);
        end
        
        keystream_chunk = typecast(state(:), 'uint8');
        
        % Calculate input/output indices
        idx_start = (b_start - 1) * BLOCK_SIZE + 1;
        idx_end = min(idx_start + length(keystream_chunk) - 1, num_bytes);
        
        % Truncate keystream if at the end of data
        len_needed = idx_end - idx_start + 1;
        keystream_chunk = keystream_chunk(1:len_needed);
        
        output(idx_start:idx_end) = bitxor(input(idx_start:idx_end), keystream_chunk);
    end
end

function x = quarter_round_vect(x, a, b, c, d)
% Optimized Vectorized ARX with MODULO FIX
% All additions must use: uint32(mod(double(A) + double(B), 4294967296))

    % 1. a += b; d ^= a; d <<<= 16
    x(a, :) = uint32(mod(double(x(a, :)) + double(x(b, :)), 4294967296));
    x(d, :) = bitxor(x(d, :), x(a, :));
    val = x(d, :);
    x(d, :) = bitor(bitshift(val, 16), bitshift(val, -16));

    % 2. c += d; b ^= c; b <<<= 12
    x(c, :) = uint32(mod(double(x(c, :)) + double(x(d, :)), 4294967296));
    x(b, :) = bitxor(x(b, :), x(c, :));
    val = x(b, :);
    x(b, :) = bitor(bitshift(val, 12), bitshift(val, -20));

    % 3. a += b; d ^= a; d <<<= 8
    x(a, :) = uint32(mod(double(x(a, :)) + double(x(b, :)), 4294967296));
    x(d, :) = bitxor(x(d, :), x(a, :));
    val = x(d, :);
    x(d, :) = bitor(bitshift(val, 8), bitshift(val, -24));

    % 4. c += d; b ^= c; b <<<= 7
    x(c, :) = uint32(mod(double(x(c, :)) + double(x(d, :)), 4294967296));
    x(b, :) = bitxor(x(b, :), x(c, :));
    val = x(b, :);
    x(b, :) = bitor(bitshift(val, 7), bitshift(val, -25));
end

function u32 = typecast_bytes_to_uint32(bytes)
    [~, ~, endian] = computer;
    u32 = typecast(bytes, 'uint32');
    if endian == 'B'
        u32 = swapbytes(u32);
    end
end
