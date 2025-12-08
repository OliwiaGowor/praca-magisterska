%% FullImageEncryption.m
% Implements the Secure Image Encryption Scheme based on Hyperchaotic System
% and 2D Compressed Sensing (Algorithms 3 & 4)
% Source: Entropy 2024, 26, 603

clc; clear; close all;

%% 1. Setup and Initialization
% Load a sample image (grayscale)
P = imread('cameraman.tif'); 
if size(P, 3) == 3
    P = rgb2gray(P);
end
P = double(P); % Work in double for calculations
[M, N] = size(P);

disp('--------------------------------------------------');
disp('Encryption Scheme Started');
disp(['Image Size: ', num2str(M), 'x', num2str(N)]);

% System Parameters (from Section 2.7 / 4.1) [cite: 357, 483]
params.a = 20;
params.b = 30;
params.x0 = 0.123987;
params.y0 = 0.987321;
params.k = 123.456;      % Intermediate key (double)
params.Pre = uint8(50);  % Intermediate key (8-bit integer)

%% 2. Encryption Process (Algorithm 3)
disp('Encrypting...');
tic;
[C, KeySeq_X, KeySeq_Y] = EncryptImage(P, params);
encryptionTime = toc;
disp(['Encryption finished in ', num2str(encryptionTime), ' seconds.']);

%% 3. Decryption Process (Algorithm 4)
disp('Decrypting...');
tic;
P_recovered = DecryptImage(C, KeySeq_X, KeySeq_Y, params);
decryptionTime = toc;
disp(['Decryption finished in ', num2str(decryptionTime), ' seconds.']);

%% 4. Verification and Visualization
% Check if perfectly reconstructed
diff = sum(sum(abs(P - P_recovered)));
if diff == 0
    disp('SUCCESS: Decrypted image matches original perfectly.');
else
    disp(['WARNING: Difference detected (Error Sum: ', num2str(diff), ')']);
end

% Display Results
figure('Name', 'Encryption/Decryption Results', 'NumberTitle', 'off');
subplot(1, 3, 1);
imshow(uint8(P));
title('Original Plaintext');

subplot(1, 3, 2);
imshow(uint8(C));
title('Encrypted Ciphertext');

subplot(1, 3, 3);
imshow(uint8(P_recovered));
title('Decrypted Image');


%% ---------------------------------------------------------
%  FUNCTION DEFINITIONS
%  ---------------------------------------------------------

% --- 2D Hyperchaotic Map (Eq. 1)  ---
function [x_next, y_next] = HyperchaoticMap(x, y, a, b)
    % x_n+1 = sin^2(a*pi*x_n + b*y_n)
    x_next = sin(a * pi * x + b * y)^2;
    
    % y_n+1 = cos^2(b*pi/y_n + a*x_n)
    % Handle potential division by zero for stability (though rare in floats)
    if y == 0
        term = 0; 
    else
        term = b * pi / y;
    end
    y_next = cos(term + a * x)^2;
end

% --- Algorithm 3: Encryption  ---
function [C, KeySeq_X, KeySeq_Y] = EncryptImage(P, params)
    [M, N] = size(P);
    
    % Initialize outputs
    C = zeros(M, N);
    KeySeq_X = zeros(1, M*N); % Store valid x keys
    KeySeq_Y = zeros(1, M*N); % Store valid y keys
    
    % Initialize State
    x = params.x0;
    y = params.y0;
    Pre = params.Pre; % 8-bit unsigned integer
    k = params.k;     % Double precision
    a = params.a;
    b = params.b;
    
    % Flag matrix to track filled positions [cite: 492]
    flag = zeros(M, N);
    
    count = 0;
    
    % Loop 1:M, 1:N [cite: 495-496]
    for i = 1:M
        for j = 1:N
            count = count + 1;
            
            % --- Find a unique empty position ---
            % "while flag(i_next, j_next) == 1 do Repeat" [cite: 503-505]
            valid_position = false;
            
            while ~valid_position
                [x, y] = HyperchaoticMap(x, y, a, b); % Iterative update [cite: 500]
                
                % Calculate coordinates (1-based index for MATLAB)
                i_next = mod(floor(x * 10^6), M) + 1; % [cite: 501]
                j_next = mod(floor(y * 10^6), N) + 1; % [cite: 502]
                
                if flag(i_next, j_next) == 0
                    valid_position = true;
                end
            end
            
            % --- Encryption Equation [cite: 520] ---
            % C(next) = bitxor(mod(P(i,j) + k, 256), Pre)
            % Note: P is double, Pre is uint8. 
            % We must be careful with types for bitxor.
            
            val_shifted = mod(P(i, j) + k, 256);
            
            % Perform XOR (Convert to uint8 for bitwise op)
            c_val = bitxor(uint8(val_shifted), Pre);
            
            % Assign to Ciphertext matrix
            C(i_next, j_next) = double(c_val);
            
            % --- Update State [cite: 521-522] ---
            Pre = c_val;           % Update Pre for next pixel diffusion
            flag(i_next, j_next) = 1; % Mark position as filled
            
            % --- Record Keys [cite: 523] ---
            KeySeq_X(count) = x;
            KeySeq_Y(count) = y;
        end
    end
end

% --- Algorithm 4: Decryption  ---
function P_recovered = DecryptImage(C, KeySeq_X, KeySeq_Y, params)
    [M, N] = size(C);
    P_recovered = zeros(M, N);
    
    % Initialize State
    Pre = params.Pre; % Must match initial encryption Pre [cite: 530]
    k = params.k;
    
    num = 1; % Key sequence counter [cite: 533]
    
    % Loop 1:M, 1:N [cite: 534-535]
    for i = 1:M
        for j = 1:N
            % --- Retrieve Coordinates from Keys [cite: 538-539] ---
            % We use the EXACT keys recorded during encryption to hit the
            % same coordinates in the same order.
            x = KeySeq_X(num);
            y = KeySeq_Y(num);
            
            i_next = mod(floor(x * 10^6), M) + 1;
            j_next = mod(floor(y * 10^6), N) + 1;
            
            % --- Decryption Equation [cite: 540] ---
            % P(i,j) = mod(bitxor(C(next), Pre) - k, 256)
            
            c_val = uint8(C(i_next, j_next)); % Retrieve ciphertext pixel
            
            % XOR with previous ciphertext (Pre) to reverse diffusion
            xor_result = double(bitxor(c_val, Pre));
            
            % Subtract k and handle modulo 256
            p_val = mod(xor_result - k, 256);
            
            P_recovered(i, j) = p_val;
            
            % --- Update State [cite: 550] ---
            Pre = c_val; % Update Pre using the CURRENT ciphertext pixel
            num = num + 1; % [cite: 551]
        end
    end
end
