function ciphertext = aes_cipher_vectorized(plaintext, w, Nr, s_box)
% AES_CIPHER_VECTORIZED - Wektoryzowane szyfrowanie AES (ECB)
% Przetwarza wszystkie bloki jednocześnie.

    % 1. Przygotowanie danych (Reshape do 4x4xN)
    % plaintext musi być wektorem uint8 o długości wielokrotności 16
    num_blocks = length(plaintext) / 16;
    state = reshape(uint8(plaintext), 4, 4, num_blocks);
    Nb = 4;

    % --- Round 0 (AddRoundKey) ---
    round_key = w(:, 1:4);
    % Implicit expansion (MATLAB R2016b+): (4x4xN) xor (4x4) działa automatycznie
    state = bitxor(state, round_key);

    % --- Rounds 1 to Nr-1 ---
    for round = 1 : Nr-1
        % 1. SubBytes (Działa na całej macierzy)
        state = s_box(double(state) + 1);
        
        % 2. ShiftRows (Operuje na wymiarze 2: kolumnach wewnątrz bloków)
        % Wiersz 1 bez zmian
        state(2, :, :) = circshift(state(2, :, :), -1, 2); % Shift left 1
        state(3, :, :) = circshift(state(3, :, :), -2, 2); % Shift left 2
        state(4, :, :) = circshift(state(4, :, :), -3, 2); % Shift left 3
        
        % 3. MixColumns (Wektoryzowane obliczenia Galois)
        % T = state * 2
        t = bitshift(state, 1);
        % Jeśli MSB był 1, XOR z 0x1B (27)
        % idx to maska logiczna dla całej macierzy 3D
        idx = state >= 128;
        t(idx) = bitxor(t(idx), 27);
        
        t2 = t;              % x2
        t3 = bitxor(t2, state); % x3 = x2 ^ x1
        
        % Macierz MixColumns:
        % 2 3 1 1
        % 1 2 3 1
        % 1 1 2 3
        % 3 1 1 2
        
        % Tworzymy nową macierz stanu (prealokacja nie jest konieczna przy przypisaniu całego wymiaru)
        state_new = state; % Kopia dla 1 i 1
        
        % Row 1: 2*s0 ^ 3*s1 ^ 1*s2 ^ 1*s3
        state_new(1,:,:) = bitxor(bitxor(bitxor(t2(1,:,:), t3(2,:,:)), state(3,:,:)), state(4,:,:));
        % Row 2: 1*s0 ^ 2*s1 ^ 3*s2 ^ 1*s3
        state_new(2,:,:) = bitxor(bitxor(bitxor(state(1,:,:), t2(2,:,:)), t3(3,:,:)), state(4,:,:));
        % Row 3: 1*s0 ^ 1*s1 ^ 2*s2 ^ 3*s3
        state_new(3,:,:) = bitxor(bitxor(bitxor(state(1,:,:), state(2,:,:)), t2(3,:,:)), t3(4,:,:));
        % Row 4: 3*s0 ^ 1*s1 ^ 1*s2 ^ 2*s3
        state_new(4,:,:) = bitxor(bitxor(bitxor(t3(1,:,:), state(2,:,:)), state(3,:,:)), t2(4,:,:));
        
        state = state_new;

        % 4. AddRoundKey
        col_start = (round * Nb) + 1;
        round_key = w(:, col_start : col_start+3);
        state = bitxor(state, round_key);
    end

    % --- Final Round ---
    state = s_box(double(state) + 1);
    
    state(2, :, :) = circshift(state(2, :, :), -1, 2);
    state(3, :, :) = circshift(state(3, :, :), -2, 2);
    state(4, :, :) = circshift(state(4, :, :), -3, 2);
    
    col_start = (Nr * Nb) + 1;
    round_key = w(:, col_start : col_start+3);
    state = bitxor(state, round_key);

    % Flatten result (4x4xN -> 16N vector)
    ciphertext = state(:)';
end