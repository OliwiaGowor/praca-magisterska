function [out1, out2] = chaos_bit_shuffling(mode, varargin)
    % DYSPOZYTOR: Pozwala wywołać funkcje lokalne z zewnątrz
    if strcmp(mode, 'encrypt')
        [out1, out2] = encrypt_image(varargin{1});
    elseif strcmp(mode, 'decrypt')
        out1 = decrypt_image(varargin{1}, varargin{2});
        out2 = []; % decrypt zwraca tylko obraz
    else
        error('Tryb musi byc "encrypt" lub "decrypt"');
    end
end

% EXTRACTED ALGORITHM FROM: encryption_GUI.m
% Original Authors: Lazaros Moysis, Ioannis Kafetzis, et al.
% Citation: Moysis, L., Kafetzis, I., Tutueva, A., Butusov, D., & Volos, C. (2022).
% Chaos-Based Image Encryption Based on Bit Level Cubic Shuffling.

function [Encryptedim, keys] = encrypt_image(input_image)
    % --- Preprocessing ---
    % Convert to gray if RGB
    if size(input_image, 3) == 3
        A = rgb2gray(input_image);
    else
        A = input_image;
    end
    
    A = double(A);
    [rows, cols] = size(A);
    
    % Bit decomposition
    Abits = zeros(rows, cols, 8);
    for i = 1:rows
        for j = 1:cols
            Abits(i,j,1:8) = bitget(A(i,j), 1:8);
        end
    end
    
    shuffled = Abits; 
    
    % --- Key Generation / Initialization ---
    AVG = mod(sum(sum(sum(A)))/(rows*cols), 1);
    
    x(1) = AVG;
    mu1 = 900 + AVG;
    a1 = 2*pi*AVG;
    
    y(1) = mod(10^6*AVG, 1);
    mu2 = 901 + y(1);
    a2 = 2*pi*y(1);
    
    z(1) = mod(10^9*AVG, 1);
    mu3 = 902 + z(1);
    a3 = 2*pi*z(1);
    
    v(1) = cos(AVG);
    mu4 = 905 + v(1);
    a4 = 2*pi*v(1);
    
    if AVG == 0
        AVG = mod(log(rows*cols + sum(sum(sum(A)))/(rows*cols)), 1);
        x(1) = AVG;
        mu1 = 900 + AVG;
        a1 = 2*pi*AVG;
        y(1) = mod(10^6*AVG, 1);
        mu2 = 901 + y(1);
        a2 = 2*pi*y(1);
        z(1) = mod(10^9*AVG, 1);
        mu3 = 902 + z(1);
        a3 = 2*pi*z(1);
        v(1) = cos(AVG);
        mu4 = 905 + v(1);
        a4 = 2*pi*v(1);
    end
    
    % The 12 generated keys
    keys = [x(1), mu1, a1, y(1), mu2, a2, z(1), mu3, a3, v(1), mu4, a4];
    
    % --- Shuffle Rows ---
    for R = 1:rows
        if R ~= 1
            x(1) = cos(mu1*(x(i)^3+x(i))+a1); 
        end
        i = 1;
        j = 0;
        per = [];
        per = floor(cols*abs(x(1)))+1;
        
        while j < cols
            i = i + 1;
            x(i) = cos(mu1*(x(i-1)^3+x(i-1))+a1);
            pos = floor(cols*abs(x(i)))+1;
            
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        left = zeros(cols, cols);
        for cc = 1:cols
            left(cc, per(cc)) = 1;
        end
        
        % Logic for the right matrix (8x8)
        clear pos
        j = 0;
        per = [];
        per = floor(8*abs(x(1)))+1;
        while j < 8
            i = i + 1;
            x(i) = cos(mu1*(x(i-1)^3+x(i-1))+a1);
            pos = floor(8*abs(x(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        right = zeros(8, 8);
        for cc = 1:8
            right(per(cc), cc) = 1;
        end
        
        shuffled(R,:,:) = left * squeeze(shuffled(R,:,:)) * right;
    end
    
    % --- Shuffle Columns ---
    for C = 1:cols
        if C ~= 1
            y(1) = cos(mu2*(y(i)^3+y(i))+a2);
        end
        i = 1;
        j = 0;
        per = [];
        per = floor(rows*abs(y(1)))+1;
        
        while j < rows
            i = i + 1;
            y(i) = cos(mu2*(y(i-1)^3+y(i-1))+a2);
            pos = floor(rows*abs(y(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        left = zeros(rows, rows);
        for cc = 1:rows
            left(cc, per(cc)) = 1;
        end
        
        clear pos
        j = 0;
        per = [];
        per = floor(8*abs(y(1)))+1;
        while j < 8
            i = i + 1;
            y(i) = cos(mu2*(y(i-1)^3+y(i-1))+a2);
            pos = floor(8*abs(y(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        right = zeros(8, 8);
        for cc = 1:8
            right(per(cc), cc) = 1;
        end
        
        shuffled(:,C,:) = left * squeeze(shuffled(:,C,:)) * right;
    end
    
    % --- Shuffle Pixel Level ---
    for P = 1:8
        if P ~= 1
            z(1) = cos(mu3*(z(i)^3+z(i))+a3);
        end
        i = 1;
        j = 0;
        per = [];
        per = floor(rows*abs(z(1)))+1;
        
        while j < rows
            i = i + 1;
            z(i) = cos(mu3*(z(i-1)^3+z(i-1))+a3);
            pos = floor(rows*abs(z(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        left = zeros(rows, rows);
        for cc = 1:rows
            left(cc, per(cc)) = 1;
        end
        
        clear pos
        j = 0;
        per = [];
        per = floor(cols*abs(z(1)))+1;
        while j < cols
            i = i + 1;
            z(i) = cos(mu3*(z(i-1)^3+z(i-1))+a3);
            pos = floor(cols*abs(z(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        right = zeros(cols, cols);
        for cc = 1:cols
            right(per(cc), cc) = 1;
        end
        
        shuffled(:,:,P) = left * squeeze(shuffled(:,:,P)) * right;
    end
    
    % --- Diffusion (XOR) ---
    shuffledstream = reshape(shuffled, 1, []);
    stream(1) = 1;
    for i = 2:rows*cols*8
        v(i) = cos(mu4*(v(i-1)^3+v(i-1))+a4);
        stream(i) = floor(mod(10^12*abs(v(i)+v(i-1)), 2));
    end
    
    encryptedstream = xor(shuffledstream, stream);
    
    % --- Reconstruction ---
    Encryptedim = zeros(rows, cols);
    test = reshape(encryptedstream, rows, cols, 8);
    for i = 1:rows
        for j = 1:cols
            Encryptedim(i,j) = bi2de(reshape(test(i,j,:), 1, []));
        end
    end
end

function Reconstructedim = decrypt_image(Encryptedim, keys)
    Encryptedim = double(Encryptedim);
    [rows, cols] = size(Encryptedim);
    
    % Bit decomposition
    Encryptedimbits = zeros(rows, cols, 8);
    for i = 1:rows
        for j = 1:cols
            Encryptedimbits(i,j,1:8) = bitget(Encryptedim(i,j), 1:8);
        end
    end
    encryptedstream = reshape(Encryptedimbits, 1, []);
    
    % Parse Keys
    x(1) = keys(1);
    mu1 = keys(2);
    a1 = keys(3);
    y(1) = keys(4);
    mu2 = keys(5);
    a2 = keys(6);
    z(1) = keys(7);
    mu3 = keys(8);
    a3 = keys(9);
    v(1) = keys(10);
    mu4 = keys(11);
    a4 = keys(12);
    
    % --- Reverse Diffusion (XOR) ---
    stream(1) = 1;
    for i = 2:rows*cols*8
        v(i) = cos(mu4*(v(i-1)^3+v(i-1))+a4);
        stream(i) = floor(mod(10^12*abs(v(i)+v(i-1)), 2));
    end
    
    Dencryptedstream = xor(encryptedstream, stream);
    Dshuffledim = reshape(Dencryptedstream, rows, cols, 8);
    
    % --- Reverse Shuffle: Pixel Level ---
    for P = 1:8
        if P ~= 1
            z(1) = cos(mu3*(z(i)^3+z(i))+a3); 
        end
        i = 1;
        j = 0;
        per = [];
        per = floor(rows*abs(z(1)))+1;
        
        while j < rows
            i = i + 1;
            z(i) = cos(mu3*(z(i-1)^3+z(i-1))+a3);
            pos = floor(rows*abs(z(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        left = zeros(rows, rows);
        for cc = 1:rows
            left(cc, per(cc)) = 1; 
        end
        
        clear pos
        j = 0;
        per = [];
        per = floor(cols*abs(z(1)))+1;
        while j < cols
            i = i + 1;
            z(i) = cos(mu3*(z(i-1)^3+z(i-1))+a3);
            pos = floor(cols*abs(z(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        right = zeros(cols, cols);
        for cc = 1:cols
            right(per(cc), cc) = 1;
        end
        
        % Using inv() as per original algorithm
        Dshuffledim(:,:,P) = inv(left) * squeeze(Dshuffledim(:,:,P)) * inv(right);
    end
    
    % --- Reverse Shuffle: Columns ---
    for C = 1:cols
        if C ~= 1
            y(1) = cos(mu2*(y(i)^3+y(i))+a2);
        end
        i = 1;
        j = 0;
        per = [];
        per = floor(rows*abs(y(1)))+1;
        
        while j < rows
            i = i + 1;
            y(i) = cos(mu2*(y(i-1)^3+y(i-1))+a2);
            pos = floor(rows*abs(y(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        left = zeros(rows, rows);
        for cc = 1:rows
            left(cc, per(cc)) = 1;
        end
        
        clear pos
        j = 0;
        per = [];
        per = floor(8*abs(y(1)))+1;
        while j < 8
            i = i + 1;
            y(i) = cos(mu2*(y(i-1)^3+y(i-1))+a2);
            pos = floor(8*abs(y(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        right = zeros(8, 8);
        for cc = 1:8
            right(per(cc), cc) = 1;
        end
        
        Dshuffledim(:,C,:) = inv(left) * squeeze(Dshuffledim(:,C,:)) * inv(right);
    end
    
    % --- Reverse Shuffle: Rows ---
    for R = 1:rows
        if R ~= 1
            x(1) = cos(mu1*(x(i)^3+x(i))+a1);
        end
        i = 1;
        j = 0;
        per = [];
        per = floor(cols*abs(x(1)))+1;
        
        while j < cols
            i = i + 1;
            x(i) = cos(mu1*(x(i-1)^3+x(i-1))+a1);
            pos = floor(cols*abs(x(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        left = zeros(cols, cols);
        for cc = 1:cols
            left(cc, per(cc)) = 1;
        end
        
        clear pos
        j = 0;
        per = [];
        per = floor(8*abs(x(1)))+1;
        while j < 8
            i = i + 1;
            x(i) = cos(mu1*(x(i-1)^3+x(i-1))+a1);
            pos = floor(8*abs(x(i)))+1;
            if per ~= pos
                j = j + 1;
                per(j) = pos;
            end
        end
        
        right = zeros(8, 8);
        for cc = 1:8
            right(per(cc), cc) = 1;
        end
        
        Dshuffledim(R,:,:) = inv(left) * squeeze(Dshuffledim(R,:,:)) * inv(right);
    end
    
    % --- Reconstruction ---
    Reconstructedim = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            Reconstructedim(i,j) = bi2de(reshape(Dshuffledim(i,j,:), 1, []));
        end
    end
    
    Reconstructedim = uint8(Reconstructedim);
end


% =========================================================
% FUNKCJA POMOCNICZA (Zastępuje Communications Toolbox)
% =========================================================
function d = bi2de_custom(b)
    % Konwertuje wektor bitów na liczbę dziesiętną.
    % Zakłada, że bity są w kolejności [LSB ... MSB] (tak jak w bitget)
    
    b = double(b);
    powers = 2 .^ (0 : length(b)-1);
    d = sum(b .* powers);
end
