function ciphertext = aes_cipher_core(plaintext, w, Nr, s_box)
% AES_CIPHER_CORE Fast AES block encryption using precomputed schedule
%
%   ciphertext = aes_cipher_core(plaintext, w, Nr, s_box)

    Nb = 4;
    state = reshape(uint8(plaintext), 4, 4);

    % --- Round 0 ---
    state = bitxor(state, w(:, 1:4));

    % --- Rounds 1 to Nr-1 ---
    % Vectorized MixColumns and Lookups
    for round = 1 : Nr-1
        % 1. SubBytes (Matrix indexing)
        state = s_box(double(state) + 1);
        
        % 2. ShiftRows 
        % Row 1 (Index 2) shift left 1
        state(2, [1 2 3 4]) = state(2, [2 3 4 1]);
        % Row 2 (Index 3) shift left 2
        state(3, [1 2 3 4]) = state(3, [3 4 1 2]);
        % Row 3 (Index 4) shift left 3
        state(4, [1 2 3 4]) = state(4, [4 1 2 3]);
        
        % 3. MixColumns 
        % T = state * 2
        t = bitshift(state, 1);
        % If MSB was 1, XOR with 0x1B
        idx = state >= 128;
        t(idx) = bitxor(t(idx), 27);
        
        % T2 = state * 2, T3 = state * 3 (which is T2 ^ state)
        t2 = t;
        t3 = bitxor(t2, state);
        
        state_out = zeros(4, 4, 'uint8');
        state_out(1,:) = bitxor(bitxor(t2(1,:), t3(2,:)), bitxor(state(3,:), state(4,:)));
        state_out(2,:) = bitxor(bitxor(state(1,:), t2(2,:)), bitxor(t3(3,:), state(4,:)));
        state_out(3,:) = bitxor(bitxor(state(1,:), state(2,:)), bitxor(t2(3,:), t3(4,:)));
        state_out(4,:) = bitxor(bitxor(t3(1,:), state(2,:)), bitxor(state(3,:), t2(4,:)));
        state = state_out;

        % 4. AddRoundKey
        col_start = (round * Nb) + 1;
        state = bitxor(state, w(:, col_start : col_start+3));
    end

    % --- Final Round ---
    state = s_box(double(state) + 1);
    
    state(2, [1 2 3 4]) = state(2, [2 3 4 1]);
    state(3, [1 2 3 4]) = state(3, [3 4 1 2]);
    state(4, [1 2 3 4]) = state(4, [4 1 2 3]);
    
    col_start = (Nr * Nb) + 1;
    state = bitxor(state, w(:, col_start : col_start+3));

    ciphertext = state(:)';
end
