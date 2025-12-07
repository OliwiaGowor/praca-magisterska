function plaintext = aes_inv_cipher_vectorized(ciphertext, w, Nr, inv_s_box, m9, m11, m13, m14)
% AES_INV_CIPHER_VECTORIZED Decrypts N blocks simultaneously
%
%   Input:
%       ciphertext: Vector of length 16*N
%       w, Nr, inv_s_box: Standard schedules
%       m9..m14: Lookup tables from get_mix_col_tables (optional, can be internal)

    % Determine number of blocks
    len = length(ciphertext);
    num_blocks = len / 16;
    Nb = 4;
    
    % Reshape to 4 x 4 x N (3D Matrix for efficient ShiftRows)
    state = reshape(uint8(ciphertext), 4, 4, num_blocks);

    % --- Initial Round (Round Nr) ---
    col_start = (Nr * Nb) + 1;
    round_key = w(:, col_start : col_start+3);
    
    % Vectorized AddRoundKey (Implicit Expansion adds 4x4 key to 4x4xN state)
    state = bitxor(state, round_key);

    % --- Rounds Nr-1 down to 1 ---
    for round = (Nr - 1) : -1 : 1
        % 1. InvShiftRows (Operates on 2nd Dimension of 3D array)
        state(2, :, :) = circshift(state(2, :, :), [0, 1]); % Shift Right 1
        state(3, :, :) = circshift(state(3, :, :), [0, 2]); % Shift Right 2
        state(4, :, :) = circshift(state(4, :, :), [0, 3]); % Shift Right 3

        % 2. InvSubBytes
        state = inv_s_box(double(state) + 1);

        % 3. AddRoundKey
        col_start = (round * Nb) + 1;
        round_key = w(:, col_start : col_start+3);
        state = bitxor(state, round_key);

        % 4. InvMixColumns (Using Lookup Tables)
        % Reshape to 2D (4 x 4N) for matrix logic
        state_2d = reshape(state, 4, []);
        
        % Lookup values (Vectorized Indexing)
        % Note: Casting to double for indexing is required
        s0 = state_2d(1,:); s1 = state_2d(2,:); s2 = state_2d(3,:); s3 = state_2d(4,:);
        
        % Table lookups replace the math loops
        v0 = double(s0) + 1; v1 = double(s1) + 1; v2 = double(s2) + 1; v3 = double(s3) + 1;
        
        % Row 0 = 0e*s0 ^ 0b*s1 ^ 0d*s2 ^ 09*s3
        state_2d(1,:) = bitxor(bitxor(bitxor(m14(v0), m11(v1)), m13(v2)), m9(v3));
        state_2d(2,:) = bitxor(bitxor(bitxor(m9(v0),  m14(v1)), m11(v2)), m13(v3));
        state_2d(3,:) = bitxor(bitxor(bitxor(m13(v0), m9(v1)),  m14(v2)), m11(v3));
        state_2d(4,:) = bitxor(bitxor(bitxor(m11(v0), m13(v1)), m9(v2)),  m14(v3));
        
        % Reshape back to 3D for next ShiftRows
        state = reshape(state_2d, 4, 4, num_blocks);
    end

    % --- Final Round ---
    state(2, :, :) = circshift(state(2, :, :), [0, 1]);
    state(3, :, :) = circshift(state(3, :, :), [0, 2]);
    state(4, :, :) = circshift(state(4, :, :), [0, 3]);

    state = inv_s_box(double(state) + 1);

    round_key = w(:, 1:4);
    state = bitxor(state, round_key);

    % Flatten result
    plaintext = state(:)';
end