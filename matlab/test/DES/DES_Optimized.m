classdef DES_Optimized < handle
    % DES_Optimized - Optimized Version (SP-Box + Permutation Lookups)

    properties (Constant)
        PC1 = [57 49 41 33 25 17 9 1 58 50 42 34 26 18 10 2 59 51 43 35 27 19 11 3 60 52 44 36 63 55 47 39 31 23 15 7 62 54 46 38 30 22 14 6 61 53 45 37 29 21 13 5 28 20 12 4];
        PC2 = [14 17 11 24 1 5 3 28 15 6 21 10 23 19 12 4 26 8 16 7 27 20 13 2 41 52 31 37 47 55 30 40 51 45 33 48 44 49 39 56 34 53 46 42 50 36 29 32];
        SHIFT_SCHEDULE = [1 1 2 2 2 2 2 2 1 2 2 2 2 2 2 1];

        IP_Table = [58 50 42 34 26 18 10 2 60 52 44 36 28 20 12 4 62 54 46 38 30 22 14 6 64 56 48 40 32 24 16 8 57 49 41 33 25 17 9 1 59 51 43 35 27 19 11 3 61 53 45 37 29 21 13 5 63 55 47 39 31 23 15 7];
        FP_Table = [40 8 48 16 56 24 64 32 39 7 47 15 55 23 63 31 38 6 46 14 54 22 62 30 37 5 45 13 53 21 61 29 36 4 44 12 52 20 60 28 35 3 43 11 51 19 59 27 34 2 42 10 50 18 58 26 33 1 41 9 49 17 57 25];

        S_RAW = [
            14 4 13 1 2 15 11 8 3 10 6 12 5 9 0 7 0 15 7 4 14 2 13 1 10 6 12 11 9 5 3 8 4 1 14 8 13 6 2 11 15 12 9 7 3 10 5 0 15 12 8 2 4 9 1 7 5 11 3 14 10 0 6 13;
            15 1 8 14 6 11 3 4 9 7 2 13 12 0 5 10 3 13 4 7 15 2 8 14 12 0 1 10 6 9 11 5 0 14 7 11 10 4 13 1 5 8 12 6 9 3 2 15 13 8 10 1 3 15 4 2 11 6 7 12 0 5 14 9;
            10 0 9 14 6 3 15 5 1 13 12 7 11 4 2 8 13 7 0 9 3 4 6 10 2 8 5 14 12 11 15 1 13 6 4 9 8 15 3 0 11 1 2 12 5 10 14 7 1 10 13 0 6 9 8 7 4 15 14 3 11 5 2 12;
            7 13 14 3 0 6 9 10 1 2 8 5 11 12 4 15 13 8 11 5 6 15 0 3 4 7 2 12 1 10 14 9 10 6 9 0 12 11 7 13 15 1 3 14 5 2 8 4 3 15 0 6 10 1 13 8 9 4 5 11 12 7 2 14;
            2 12 4 1 7 10 11 6 8 5 3 15 13 0 14 9 14 11 2 12 4 7 13 1 5 0 15 10 3 9 8 6 4 2 1 11 10 13 7 8 15 9 12 5 6 3 0 14 11 8 12 7 1 14 2 13 6 15 0 9 10 4 5 3;
            12 1 10 15 9 2 6 8 0 13 3 4 14 7 5 11 10 15 4 2 7 12 9 5 6 1 13 14 0 11 3 8 9 14 15 5 2 8 12 3 7 0 4 10 1 13 11 6 4 3 2 12 9 5 15 10 11 14 1 7 6 0 8 13;
            4 11 2 14 15 0 8 13 3 12 9 7 5 10 6 1 13 0 11 7 4 9 1 10 14 3 5 12 2 15 8 6 1 4 11 13 12 3 7 14 10 15 6 8 0 5 9 2 6 11 13 8 1 4 10 7 9 5 0 15 14 2 3 12;
            13 2 8 4 6 15 11 1 10 9 3 14 5 0 12 7 1 15 13 8 10 3 7 4 12 5 6 11 0 14 9 2 7 11 4 1 9 12 14 2 0 6 10 13 15 3 5 8 2 1 14 7 4 10 8 13 15 12 9 0 3 5 6 11
            ];
        P_PERM = [16 7 20 21 29 12 28 17 1 15 23 26 5 18 31 10 2 8 24 14 32 27 3 9 19 13 30 6 22 11 4 25];
    end

    properties
        SP_BOXES % 8 x 64 uint32
        IP_LUT   % 8 x 256 uint64 (Fast Initial Permutation)
        FP_LUT   % 8 x 256 uint64 (Fast Final Permutation)
    end

    methods
        function obj = DES_Optimized()
            obj.init_sp_tables();
            obj.init_perm_lookups();
        end

        function init_sp_tables(obj)
            % Precompute S-Box + P Permutation
            obj.SP_BOXES = zeros(8, 64, 'uint32');
            for k = 1:8
                for x = 0:63
                    b1 = bitget(x, 6);
                    b6 = bitget(x, 1);
                    row = b1*2 + b6 + 1;
                    col = bitand(bitshift(x, -1), 15) + 1;
                    s_val = obj.S_RAW(k, (row-1)*16 + col);

                    val_P = uint32(0);
                    for bit_idx = 1:4
                        if bitget(s_val, 5 - bit_idx)
                            input_pos_P = (k-1)*4 + bit_idx;
                            target_pos = find(obj.P_PERM == input_pos_P);
                            val_P = bitset(val_P, 33 - target_pos);
                        end
                    end
                    obj.SP_BOXES(k, x+1) = val_P;
                end
            end
        end

        function init_perm_lookups(obj)
            % Creates LUT tables for IP and FP.
            % Instead of permuting 64 bits individually, permute 8 bytes.
            % For each byte (0-255), we know exactly where its bits will end up.

            obj.IP_LUT = zeros(8, 256, 'uint64');
            obj.FP_LUT = zeros(8, 256, 'uint64');

            for byteIdx = 1:8
                % Shift to place the test byte in the correct uint64 position
                % byteIdx 1 = bits 64..57 (MSB in FIPS), byteIdx 8 = bits 8..1
                shiftAmt = (8 - byteIdx) * 8;

                for val = 0:255
                    % Simulate input with only one byte set
                    testInput = bitshift(uint64(val), shiftAmt);

                    % Calculate where bits land
                    perm_IP = obj.slow_perm(testInput, obj.IP_Table);
                    perm_FP = obj.slow_perm(testInput, obj.FP_Table);

                    obj.IP_LUT(byteIdx, val + 1) = perm_IP;
                    obj.FP_LUT(byteIdx, val + 1) = perm_FP;
                end
            end
        end

        function cipher64 = encrypt_data(obj, plain64, key64, iv64)
            subKeys = obj.generate_subkeys(key64);
            numBlocks = length(plain64);
            cipher64 = zeros(size(plain64), 'uint64');
            prevBlock = iv64;
            mask32 = uint64(4294967295);

            % --- OPTYMALIZACJA: Cache'owanie zmiennych lokalnych ---
            % Dostęp do zmiennych lokalnych jest szybszy niż do pól obiektu (obj.X)
            local_IP_LUT = obj.IP_LUT;
            local_FP_LUT = obj.FP_LUT;
            local_SP_BOXES = obj.SP_BOXES;
            % -------------------------------------------------------

            for i = 1:numBlocks
                block = bitxor(plain64(i), prevBlock);

                % IP Permutation (Lookup)
                block = bitxor(bitxor(bitxor( ...
                    local_IP_LUT(1, bitshift(block, -56) + 1),       ...
                    local_IP_LUT(2, bitand(bitshift(block, -48), 255) + 1)), ...
                    bitxor( ...
                    local_IP_LUT(3, bitand(bitshift(block, -40), 255) + 1),  ...
                    local_IP_LUT(4, bitand(bitshift(block, -32), 255) + 1))), ...
                    bitxor(bitxor( ...
                    local_IP_LUT(5, bitand(bitshift(block, -24), 255) + 1),  ...
                    local_IP_LUT(6, bitand(bitshift(block, -16), 255) + 1)), ...
                    bitxor( ...
                    local_IP_LUT(7, bitand(bitshift(block, -8),  255) + 1),  ...
                    local_IP_LUT(8, bitand(block, 255) + 1))));

                L = uint32(bitshift(block, -32));
                R = uint32(bitand(block, mask32));

                % 16 Feistel Rounds
                for r = 1:16
                    L_next = R;
                    K_curr = subKeys(r);
                    f_val = uint32(0);

                    % Użycie lokalnej kopii SP_BOXES
                    idx = bitxor(bitand(bitshift(R, -27) + bitshift(R, 5), 63), uint32(bitand(bitshift(K_curr, -42), 63)));
                    f_val = bitxor(f_val, local_SP_BOXES(1, idx + 1));
                    idx = bitxor(bitand(bitshift(R, -23), 63), uint32(bitand(bitshift(K_curr, -36), 63)));
                    f_val = bitxor(f_val, local_SP_BOXES(2, idx + 1));
                    idx = bitxor(bitand(bitshift(R, -19), 63), uint32(bitand(bitshift(K_curr, -30), 63)));
                    f_val = bitxor(f_val, local_SP_BOXES(3, idx + 1));
                    idx = bitxor(bitand(bitshift(R, -15), 63), uint32(bitand(bitshift(K_curr, -24), 63)));
                    f_val = bitxor(f_val, local_SP_BOXES(4, idx + 1));
                    idx = bitxor(bitand(bitshift(R, -11), 63), uint32(bitand(bitshift(K_curr, -18), 63)));
                    f_val = bitxor(f_val, local_SP_BOXES(5, idx + 1));
                    idx = bitxor(bitand(bitshift(R, -7), 63), uint32(bitand(bitshift(K_curr, -12), 63)));
                    f_val = bitxor(f_val, local_SP_BOXES(6, idx + 1));
                    idx = bitxor(bitand(bitshift(R, -3), 63), uint32(bitand(bitshift(K_curr, -6), 63)));
                    f_val = bitxor(f_val, local_SP_BOXES(7, idx + 1));
                    idx = bitxor(bitand(bitshift(R, 1) + bitshift(R, -31), 63), uint32(bitand(K_curr, 63)));
                    f_val = bitxor(f_val, local_SP_BOXES(8, idx + 1));

                    R = bitxor(L, f_val);
                    L = L_next;
                end

                preOut = bitor(bitshift(uint64(R), 32), uint64(L));

                % FP Permutation
                outBlock = bitxor(bitxor(bitxor( ...
                    local_FP_LUT(1, bitshift(preOut, -56) + 1),       ...
                    local_FP_LUT(2, bitand(bitshift(preOut, -48), 255) + 1)), ...
                    bitxor( ...
                    local_FP_LUT(3, bitand(bitshift(preOut, -40), 255) + 1),  ...
                    local_FP_LUT(4, bitand(bitshift(preOut, -32), 255) + 1))), ...
                    bitxor(bitxor( ...
                    local_FP_LUT(5, bitand(bitshift(preOut, -24), 255) + 1),  ...
                    local_FP_LUT(6, bitand(bitshift(preOut, -16), 255) + 1)), ...
                    bitxor( ...
                    local_FP_LUT(7, bitand(bitshift(preOut, -8),  255) + 1),  ...
                    local_FP_LUT(8, bitand(preOut, 255) + 1))));

                cipher64(i) = outBlock;
                prevBlock = outBlock;
            end
        end
        function plain64 = decrypt_data(obj, cipher64, key64, iv64)
            % --- OPTYMALIZACJA: PEŁNA WEKTORYZACJA ---
            % W CBC deszyfrowanie nie musi być sekwencyjne.
            % Możemy obliczyć funkcję DES^-1 dla wszystkich bloków jednocześnie.

            subKeys = obj.generate_subkeys(key64);
            mask32 = uint64(4294967295);

            % Pobieramy lokalne kopie dla szybkości
            local_IP_LUT = obj.IP_LUT;
            local_FP_LUT = obj.FP_LUT;
            local_SP_BOXES = obj.SP_BOXES;

            % 1. Initial Permutation (Wektorowo dla całej tablicy cipher64)
            % MATLAB automatycznie iteruje po wszystkich elementach wektora
            block = bitxor(bitxor(bitxor( ...
                local_IP_LUT(1, bitshift(cipher64, -56) + 1),       ...
                local_IP_LUT(2, bitand(bitshift(cipher64, -48), 255) + 1)), ...
                bitxor( ...
                local_IP_LUT(3, bitand(bitshift(cipher64, -40), 255) + 1),  ...
                local_IP_LUT(4, bitand(bitshift(cipher64, -32), 255) + 1))), ...
                bitxor(bitxor( ...
                local_IP_LUT(5, bitand(bitshift(cipher64, -24), 255) + 1),  ...
                local_IP_LUT(6, bitand(bitshift(cipher64, -16), 255) + 1)), ...
                bitxor( ...
                local_IP_LUT(7, bitand(bitshift(cipher64, -8),  255) + 1),  ...
                local_IP_LUT(8, bitand(cipher64, 255) + 1))));

            L = uint32(bitshift(block, -32));
            R = uint32(bitand(block, mask32));

            % 2. 16 Rund Feistela (Pętla tylko 16 razy, wewnątrz operacje na wektorach)
            for r = 16:-1:1
                L_next = R;
                K_curr = subKeys(r);
                f_val = uint32(0);

                % Obliczamy indeksy S-Box dla CAŁEGO wektora R naraz
                idx = bitxor(bitand(bitshift(R, -27) + bitshift(R, 5), 63), uint32(bitand(bitshift(K_curr, -42), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(1, idx + 1)); % Indeksowanie wektorem działa w MATLAB b. szybko

                idx = bitxor(bitand(bitshift(R, -23), 63), uint32(bitand(bitshift(K_curr, -36), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(2, idx + 1));

                idx = bitxor(bitand(bitshift(R, -19), 63), uint32(bitand(bitshift(K_curr, -30), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(3, idx + 1));

                idx = bitxor(bitand(bitshift(R, -15), 63), uint32(bitand(bitshift(K_curr, -24), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(4, idx + 1));

                idx = bitxor(bitand(bitshift(R, -11), 63), uint32(bitand(bitshift(K_curr, -18), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(5, idx + 1));

                idx = bitxor(bitand(bitshift(R, -7), 63), uint32(bitand(bitshift(K_curr, -12), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(6, idx + 1));

                idx = bitxor(bitand(bitshift(R, -3), 63), uint32(bitand(bitshift(K_curr, -6), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(7, idx + 1));

                idx = bitxor(bitand(bitshift(R, 1) + bitshift(R, -31), 63), uint32(bitand(K_curr, 63)));
                f_val = bitxor(f_val, local_SP_BOXES(8, idx + 1));

                R = bitxor(L, f_val);
                L = L_next;
            end

            preOut = bitor(bitshift(uint64(R), 32), uint64(L));

            % 3. Final Permutation (Wektorowo)
            decryptedBlock = bitxor(bitxor(bitxor( ...
                local_FP_LUT(1, bitshift(preOut, -56) + 1),       ...
                local_FP_LUT(2, bitand(bitshift(preOut, -48), 255) + 1)), ...
                bitxor( ...
                local_FP_LUT(3, bitand(bitshift(preOut, -40), 255) + 1),  ...
                local_FP_LUT(4, bitand(bitshift(preOut, -32), 255) + 1))), ...
                bitxor(bitxor( ...
                local_FP_LUT(5, bitand(bitshift(preOut, -24), 255) + 1),  ...
                local_FP_LUT(6, bitand(bitshift(preOut, -16), 255) + 1)), ...
                bitxor( ...
                local_FP_LUT(7, bitand(bitshift(preOut, -8),  255) + 1),  ...
                local_FP_LUT(8, bitand(preOut, 255) + 1))));

            % 4. CBC XOR Recovery
            % Plain[i] = Dec(Cipher[i]) XOR Cipher[i-1]
            % Wektoryzacja: Przesuwamy wektor cipher64 i doklejamy IV na początek
            prevBlocks = [iv64, cipher64(1:end-1)];
            plain64 = bitxor(decryptedBlock, prevBlocks);
        end

        function subKeys = generate_subkeys(obj, key64)
            % It runs only once per image
            K56 = obj.slow_perm_56(key64, obj.PC1);
            C = uint32(bitshift(K56, -28));
            D = uint32(bitand(K56, 268435455));
            subKeys = zeros(16, 1, 'uint64');
            for i = 1:16
                shift = obj.SHIFT_SCHEDULE(i);
                C = bitand(bitor(bitshift(C, shift), bitshift(C, shift-28)), 268435455);
                D = bitand(bitor(bitshift(D, shift), bitshift(D, shift-28)), 268435455);
                CD = bitor(bitshift(uint64(C), 28), uint64(D));
                subKeys(i) = obj.slow_perm_48(CD, obj.PC2);
            end
        end

        function out = slow_perm(obj, in, table)
            out = uint64(0);
            for i = 1:64
                if bitget(in, 65 - table(i))
                    out = bitset(out, 65 - i);
                end
            end
        end
        function out = slow_perm_56(obj, in, table)
            out = uint64(0);
            for i = 1:56
                if bitget(in, 65 - table(i))
                    out = bitset(out, 57 - i);
                end
            end
        end
        function out = slow_perm_48(obj, in, table)
            out = uint64(0);
            for i = 1:48
                if bitget(in, 57 - table(i))
                    out = bitset(out, 49 - i);
                end
            end
        end

        % Wklej to wewnątrz bloku methods w DES_Optimized.m
        function cipher64 = encrypt_data_ecb(obj, plain64, key64)
            subKeys = obj.generate_subkeys(key64);
            mask32 = uint64(4294967295);

            % Cache'owanie zmiennych dla wydajności
            local_IP_LUT = obj.IP_LUT;
            local_FP_LUT = obj.FP_LUT;
            local_SP_BOXES = obj.SP_BOXES;

            % W trybie ECB (Pikselowym) NIE XORujemy z prevBlock/IV.
            % Każdy blok plain64 (reprezentujący 8 pikseli) wchodzi bezpośrednio.
            block = plain64;

            % 1. Initial Permutation (Wektorowo)
            block = bitxor(bitxor(bitxor( ...
                local_IP_LUT(1, bitshift(block, -56) + 1),       ...
                local_IP_LUT(2, bitand(bitshift(block, -48), 255) + 1)), ...
                bitxor( ...
                local_IP_LUT(3, bitand(bitshift(block, -40), 255) + 1),  ...
                local_IP_LUT(4, bitand(bitshift(block, -32), 255) + 1))), ...
                bitxor(bitxor( ...
                local_IP_LUT(5, bitand(bitshift(block, -24), 255) + 1),  ...
                local_IP_LUT(6, bitand(bitshift(block, -16), 255) + 1)), ...
                bitxor( ...
                local_IP_LUT(7, bitand(bitshift(block, -8),  255) + 1),  ...
                local_IP_LUT(8, bitand(block, 255) + 1))));

            L = uint32(bitshift(block, -32));
            R = uint32(bitand(block, mask32));

            % 2. 16 Rund Feistela
            for r = 1:16
                L_next = R;
                K_curr = subKeys(r);
                f_val = uint32(0);

                idx = bitxor(bitand(bitshift(R, -27) + bitshift(R, 5), 63), uint32(bitand(bitshift(K_curr, -42), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(1, idx + 1));
                idx = bitxor(bitand(bitshift(R, -23), 63), uint32(bitand(bitshift(K_curr, -36), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(2, idx + 1));
                idx = bitxor(bitand(bitshift(R, -19), 63), uint32(bitand(bitshift(K_curr, -30), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(3, idx + 1));
                idx = bitxor(bitand(bitshift(R, -15), 63), uint32(bitand(bitshift(K_curr, -24), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(4, idx + 1));
                idx = bitxor(bitand(bitshift(R, -11), 63), uint32(bitand(bitshift(K_curr, -18), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(5, idx + 1));
                idx = bitxor(bitand(bitshift(R, -7), 63), uint32(bitand(bitshift(K_curr, -12), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(6, idx + 1));
                idx = bitxor(bitand(bitshift(R, -3), 63), uint32(bitand(bitshift(K_curr, -6), 63)));
                f_val = bitxor(f_val, local_SP_BOXES(7, idx + 1));
                idx = bitxor(bitand(bitshift(R, 1) + bitshift(R, -31), 63), uint32(bitand(K_curr, 63)));
                f_val = bitxor(f_val, local_SP_BOXES(8, idx + 1));

                R = bitxor(L, f_val);
                L = L_next;
            end

            preOut = bitor(bitshift(uint64(R), 32), uint64(L));

            % 3. Final Permutation
            cipher64 = bitxor(bitxor(bitxor( ...
                local_FP_LUT(1, bitshift(preOut, -56) + 1),       ...
                local_FP_LUT(2, bitand(bitshift(preOut, -48), 255) + 1)), ...
                bitxor( ...
                local_FP_LUT(3, bitand(bitshift(preOut, -40), 255) + 1),  ...
                local_FP_LUT(4, bitand(bitshift(preOut, -32), 255) + 1))), ...
                bitxor(bitxor( ...
                local_FP_LUT(5, bitand(bitshift(preOut, -24), 255) + 1),  ...
                local_FP_LUT(6, bitand(bitshift(preOut, -16), 255) + 1)), ...
                bitxor( ...
                local_FP_LUT(7, bitand(bitshift(preOut, -8),  255) + 1),  ...
                local_FP_LUT(8, bitand(preOut, 255) + 1))));
        end

        % Wklej to do bloku 'methods' w pliku testt/DES/DES_Optimized.m

function plain64 = decrypt_data_ecb(obj, cipher64, key64)
    % DECRYPT_DATA_ECB - Deszyfrowanie niezależnych bloków (tryb ECB/Pixel)
    
    subKeys = obj.generate_subkeys(key64);
    mask32 = uint64(4294967295);
    
    % Cache dla wydajności
    local_IP_LUT = obj.IP_LUT;
    local_FP_LUT = obj.FP_LUT;
    local_SP_BOXES = obj.SP_BOXES;

    % 1. Initial Permutation (Wektorowo)
    block = bitxor(bitxor(bitxor( ...
        local_IP_LUT(1, bitshift(cipher64, -56) + 1),       ...
        local_IP_LUT(2, bitand(bitshift(cipher64, -48), 255) + 1)), ...
        bitxor( ...
        local_IP_LUT(3, bitand(bitshift(cipher64, -40), 255) + 1),  ...
        local_IP_LUT(4, bitand(bitshift(cipher64, -32), 255) + 1))), ...
        bitxor(bitxor( ...
        local_IP_LUT(5, bitand(bitshift(cipher64, -24), 255) + 1),  ...
        local_IP_LUT(6, bitand(bitshift(cipher64, -16), 255) + 1)), ...
        bitxor( ...
        local_IP_LUT(7, bitand(bitshift(cipher64, -8),  255) + 1),  ...
        local_IP_LUT(8, bitand(cipher64, 255) + 1))));

    L = uint32(bitshift(block, -32));
    R = uint32(bitand(block, mask32));

    % 2. 16 Rund Feistela (Odwrotna kolejność kluczy)
    for r = 16:-1:1
        L_next = R;
        K_curr = subKeys(r);
        f_val = uint32(0);
        
        idx = bitxor(bitand(bitshift(R, -27) + bitshift(R, 5), 63), uint32(bitand(bitshift(K_curr, -42), 63)));
        f_val = bitxor(f_val, local_SP_BOXES(1, idx + 1));
        idx = bitxor(bitand(bitshift(R, -23), 63), uint32(bitand(bitshift(K_curr, -36), 63)));
        f_val = bitxor(f_val, local_SP_BOXES(2, idx + 1));
        idx = bitxor(bitand(bitshift(R, -19), 63), uint32(bitand(bitshift(K_curr, -30), 63)));
        f_val = bitxor(f_val, local_SP_BOXES(3, idx + 1));
        idx = bitxor(bitand(bitshift(R, -15), 63), uint32(bitand(bitshift(K_curr, -24), 63)));
        f_val = bitxor(f_val, local_SP_BOXES(4, idx + 1));
        idx = bitxor(bitand(bitshift(R, -11), 63), uint32(bitand(bitshift(K_curr, -18), 63)));
        f_val = bitxor(f_val, local_SP_BOXES(5, idx + 1));
        idx = bitxor(bitand(bitshift(R, -7), 63), uint32(bitand(bitshift(K_curr, -12), 63)));
        f_val = bitxor(f_val, local_SP_BOXES(6, idx + 1));
        idx = bitxor(bitand(bitshift(R, -3), 63), uint32(bitand(bitshift(K_curr, -6), 63)));
        f_val = bitxor(f_val, local_SP_BOXES(7, idx + 1));
        idx = bitxor(bitand(bitshift(R, 1) + bitshift(R, -31), 63), uint32(bitand(K_curr, 63)));
        f_val = bitxor(f_val, local_SP_BOXES(8, idx + 1));
        
        R = bitxor(L, f_val);
        L = L_next;
    end

    preOut = bitor(bitshift(uint64(R), 32), uint64(L));

    % 3. Final Permutation
    plain64 = bitxor(bitxor(bitxor( ...
        local_FP_LUT(1, bitshift(preOut, -56) + 1),       ...
        local_FP_LUT(2, bitand(bitshift(preOut, -48), 255) + 1)), ...
        bitxor( ...
        local_FP_LUT(3, bitand(bitshift(preOut, -40), 255) + 1),  ...
        local_FP_LUT(4, bitand(bitshift(preOut, -32), 255) + 1))), ...
        bitxor(bitxor( ...
        local_FP_LUT(5, bitand(bitshift(preOut, -24), 255) + 1),  ...
        local_FP_LUT(6, bitand(bitshift(preOut, -16), 255) + 1)), ...
        bitxor( ...
        local_FP_LUT(7, bitand(bitshift(preOut, -8),  255) + 1),  ...
        local_FP_LUT(8, bitand(preOut, 255) + 1))));
        
    % W trybie ECB NIE wykonujemy XOR z prevBlock!
end
    end
end
