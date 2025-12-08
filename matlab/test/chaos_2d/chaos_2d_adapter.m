function [out1, out2] = chaos_2d_adapter(mode, varargin)
    % CHAOS_2D_ADAPTER - Exact implementation of the paper's algorithm.
    % Reference: Entropy 2024, 26, 603, Algorithm 3 & 4
    
    if strcmp(mode, 'encrypt')
        [out1, out2] = encrypt_core(varargin{1}); 
    elseif strcmp(mode, 'decrypt')
        out1 = decrypt_core(varargin{1}, varargin{2}); 
        out2 = [];
    else
        error('Mode must be "encrypt" or "decrypt"');
    end
end

% =========================================================
% CORE: ENCRYPTION (Algorithm 3)
% =========================================================
function [C, keys] = encrypt_core(P)
    % 1. Parameters [cite: 357, 483, 560]
    params.a = 20;
    params.b = 30;
    params.x0 = 0.123987;
    params.y0 = 0.987321;
    params.k = 123.456;      
    params.Pre = 50;  % Kept as double for calculation

    % 2. Preprocessing
    if size(P, 3) == 3
        P = rgb2gray(P);
    end
    P = double(P);
    [M, N] = size(P);
    
    % 3. Allocation
    C = zeros(M, N);
    KeySeq_X = zeros(1, M*N); 
    KeySeq_Y = zeros(1, M*N);
    
    % 4. State
    x = params.x0;
    y = params.y0;
    Pre = double(params.Pre); % Use double for bitxor input
    k = params.k;     
    a = params.a;
    b = params.b;
    
    flag = false(M, N);
    count = 0;
    
    % 5. Encryption Loop
    for i = 1:M
        for j = 1:N
            count = count + 1;
            
            % --- Random Coordinate Generation [cite: 501-505] ---
            valid_position = false;
            while ~valid_position
                % Hyperchaotic Map [cite: 66-68]
                x_next = sin(a * pi * x + b * y)^2;
                if y == 0, term = 0; else, term = b * pi / y; end
                y_next = cos(term + a * x)^2;
                x = x_next;
                y = y_next;
                
                % Coordinates [cite: 501-502]
                i_next = mod(floor(x * 10^6), M) + 1; 
                j_next = mod(floor(y * 10^6), N) + 1; 
                
                if ~flag(i_next, j_next)
                    valid_position = true;
                end
            end
            
            % --- Encryption  ---
            % Strict formula: C = bitxor(mod(P + k, 256), Pre)
            % Use double for modulo, then cast to uint8 ONLY for bitxor
            val_shifted = mod(P(i, j) + k, 256);
            
            c_val = bitxor(uint8(val_shifted), uint8(Pre));
            
            % Store as double in C matrix
            C(i_next, j_next) = double(c_val);
            
            % --- Update State  ---
            Pre = double(c_val); % Convert back to double for next iteration storage
            flag(i_next, j_next) = true; 
            
            % Save Keys [cite: 523]
            KeySeq_X(count) = x;
            KeySeq_Y(count) = y;
        end
    end
    
    % Return uint8 for image display/saving
    C = uint8(C);
    
    keys.params = params;
    keys.KeySeq_X = KeySeq_X;
    keys.KeySeq_Y = KeySeq_Y;
end

% =========================================================
% CORE: DECRYPTION (Algorithm 4)
% =========================================================
function P_recovered = decrypt_core(C, keys)
    C = double(C);
    [M, N] = size(C);
    P_recovered = zeros(M, N);
    
    params = keys.params;
    KeySeq_X = keys.KeySeq_X;
    KeySeq_Y = keys.KeySeq_Y;
    
    Pre = double(params.Pre); 
    k = params.k;
    num = 1; 
    
    for i = 1:M
        for j = 1:N
            % Retrieve coordinates [cite: 538-539]
            x = KeySeq_X(num);
            y = KeySeq_Y(num);
            
            i_next = mod(floor(x * 10^6), M) + 1;
            j_next = mod(floor(y * 10^6), N) + 1;
            
            % --- Decryption [cite: 540] ---
            % Formula: P = mod(bitxor(C, Pre) - k, 256)
            c_val = C(i_next, j_next);
            
            xor_result = double(bitxor(uint8(c_val), uint8(Pre)));
            
            % Handle negative results in modulo
            p_val = mod(xor_result - k, 256);
            
            P_recovered(i, j) = p_val;
            
            % Update State [cite: 550]
            Pre = double(c_val); 
            num = num + 1; 
        end
    end
    
    P_recovered = uint8(P_recovered);
end
