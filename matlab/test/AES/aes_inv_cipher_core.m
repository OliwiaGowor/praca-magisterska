function plaintext = aes_inv_cipher_core(ciphertext, w, Nr, inv_s_box)
% AES_INV_CIPHER_CORE Fast AES decryption using precomputed schedule
%
%   plaintext = aes_inv_cipher_core(ciphertext, w, Nr, inv_s_box)
%
%   Input:
%       ciphertext: 1x16 uint8
%       w:          Expanded Key Schedule (from aes_init)
%       Nr:         Number of rounds (10, 12, or 14)
%       inv_s_box:  Inverse S-Box table (256 uint8s)

    Nb = 4;
    state = reshape(uint8(ciphertext), 4, 4);

    % --- Initial Round (Round Nr) ---
    % AddRoundKey using the LAST round key
    col_start = (Nr * Nb) + 1;
    round_key = w(:, col_start : col_start+3);
    state = bitxor(state, round_key);

    % --- Rounds Nr-1 down to 1 ---
    for round = (Nr - 1) : -1 : 1
        % 1. InvShiftRows
        % Row 2 shift right 1
        state(2, [1 2 3 4]) = state(2, [4 1 2 3]);
        % Row 3 shift right 2
        state(3, [1 2 3 4]) = state(3, [3 4 1 2]);
        % Row 4 shift right 3
        state(4, [1 2 3 4]) = state(4, [2 3 4 1]);

        % 2. InvSubBytes
        state = inv_s_box(double(state) + 1);

        % 3. AddRoundKey
        col_start = (round * Nb) + 1;
        round_key = w(:, col_start : col_start+3);
        state = bitxor(state, round_key);

        % 4. InvMixColumns (Vectorized)
        % We need to multiply by 0x09, 0x0b, 0x0d, 0x0e
        % 9=8+1, 11=8+2+1, 13=8+4+1, 14=8+4+2
        
        % Calculate intermediate terms x2, x4, x8
        t_1 = state;
        t_2 = aes_xtime(t_1);
        t_4 = aes_xtime(t_2);
        t_8 = aes_xtime(t_4);
        
        % Combine terms to get x9, x11, x13, x14
        t_9  = bitxor(t_8, t_1);
        t_11 = bitxor(bitxor(t_8, t_2), t_1);
        t_13 = bitxor(bitxor(t_8, t_4), t_1);
        t_14 = bitxor(bitxor(t_8, t_4), t_2);
        
        state_out = zeros(4, 4, 'uint8');
        
        % Matrix Multiply:
        % [0e 0b 0d 09]
        % [09 0e 0b 0d]
        % [0d 09 0e 0b]
        % [0b 0d 09 0e]
        
        % Row 1
        state_out(1,:) = bitxor(bitxor(bitxor(t_14(1,:), t_11(2,:)), t_13(3,:)), t_9(4,:));
        % Row 2
        state_out(2,:) = bitxor(bitxor(bitxor(t_9(1,:), t_14(2,:)), t_11(3,:)), t_13(4,:));
        % Row 3
        state_out(3,:) = bitxor(bitxor(bitxor(t_13(1,:), t_9(2,:)), t_14(3,:)), t_11(4,:));
        % Row 4
        state_out(4,:) = bitxor(bitxor(bitxor(t_11(1,:), t_13(2,:)), t_9(3,:)), t_14(4,:));
        
        state = state_out;
    end

    % --- Final Round (Round 0) ---
    % 1. InvShiftRows
    state(2, [1 2 3 4]) = state(2, [4 1 2 3]);
    state(3, [1 2 3 4]) = state(3, [3 4 1 2]);
    state(4, [1 2 3 4]) = state(4, [2 3 4 1]);

    % 2. InvSubBytes
    state = inv_s_box(double(state) + 1);

    % 3. AddRoundKey (First Key)
    round_key = w(:, 1:4);
    state = bitxor(state, round_key);

    plaintext = state(:)';
end

function out = aes_xtime(in)
    % AES 'xtime' operation (multiply by 2 in GF(2^8))
    % Vectorized for entire 4x4 matrix
    
    % Shift left
    t = bitshift(in, 1);
    
    % If MSB was 1, XOR with 0x1B (27)
    % We check the input 'in' for >= 128 (0x80)
    idx = in >= 128;
    t(idx) = bitxor(t(idx), 27);
    
    out = t;
end