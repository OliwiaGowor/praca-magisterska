% Clear workspace
clear; clc; close all;

% --- 1. Load Image ---
% Ensure you have an image named 'lena.png' or change the filename
% The paper analyzes grayscale images [cite: 161]
originalImage = imread('peppers.png'); % Using standard MATLAB demo image
if size(originalImage, 3) == 3
    originalImage = rgb2gray(originalImage);
end

imgDouble = double(originalImage);
[M, N] = size(imgDouble);

fprintf('Image Size: %dx%d\n', M, N);

% --- 2. Parameters from Paper  ---
% Two-Stage Logistic Map Parameters
x1 = 0.41;       % Initial value for map 1
x2 = 0.87;       % Initial value for map 2
lambda1 = 3.95;  % Branch parameter (gamma) for map 1
lambda2 = 3.80;  % Branch parameter (gamma) for map 2

% M-Sequence / Permutation Parameters
% The paper lists x'0 = 10, y'0 = 2. These act as seeds/shifts for the permutation.
m_seed = 10;     
iterations = 1;  % "Let the number of iterations be a" ... "gamma=1" 

% --- 3. Encryption ---
fprintf('Encrypting...\n');
encryptedImage = imgDouble;

for k = 1:iterations
    % Step A: Substitution (Diffusion) using Two-Stage Logistic Map
    encryptedImage = logistic_substitution_encrypt(encryptedImage, x1, x2, lambda1, lambda2);
    
    % Step B: Scrambling (Permutation) using M-Sequence
    encryptedImage = m_sequence_permute(encryptedImage, m_seed, 'encrypt');
end

encryptedImage = uint8(encryptedImage);

% --- 4. Decryption ---
fprintf('Decrypting...\n');
decryptedImage = double(encryptedImage);

for k = 1:iterations
    % Inverse Step B: Inverse Scrambling
    decryptedImage = m_sequence_permute(decryptedImage, m_seed, 'decrypt');
    
    % Inverse Step A: Inverse Substitution
    decryptedImage = logistic_substitution_decrypt(decryptedImage, x1, x2, lambda1, lambda2);
end

decryptedImage = uint8(decryptedImage);

% --- 5. Display Results ---
figure('Name', 'Encryption Results based on Hu & Tian (2020)', 'NumberTitle', 'off');

subplot(1, 3, 1);
imshow(originalImage);
title('Original Image');

subplot(1, 3, 2);
imshow(encryptedImage);
title('Encrypted Image');

subplot(1, 3, 3);
imshow(decryptedImage);
title('Decrypted Image');

% Verify lossless decryption
diff = sum(abs(double(originalImage(:)) - double(decryptedImage(:))));
fprintf('Total Pixel Difference: %d\n', diff);
