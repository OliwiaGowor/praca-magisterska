classdef DES_ImgProcess_Opt
    % Wrapper obsługujący konwersję obrazu przy użyciu Wektoryzacji
    % i wywołujący DES_Optimized
    
    properties (Constant)
        Engine = DES_Optimized();
    end
    
    methods (Static)
        function [cipherBytes, cipherImgVis] = encrypt(imgUint8, keyHex, ivHex)
            % 1. Wektoryzacja danych wejściowych
            % a. Padding PKCS#7 (operacje na wektorze uint8)
            flatPixels = imgUint8(:)'; 
            padLen = 8 - mod(length(flatPixels), 8);
            paddedPixels = [flatPixels, repmat(uint8(padLen), 1, padLen)];
            
            % b. Wektoryzacja: uint8 -> uint64 (Typecast)
            % MATLAB na PC jest zazwyczaj Little Endian, a DES wymaga Big Endian.
            % Używamy swapbytes do korekcji.
            blocksIn = swapbytes(typecast(paddedPixels, 'uint64'));
            
            % 2. Przygotowanie kluczy (konwersja hex -> uint64)
            key64 = DES_ImgProcess_Opt.hexToUint64(keyHex);
            iv64  = DES_ImgProcess_Opt.hexToUint64(ivHex);
            
            % 3. Szybkie Szyfrowanie (DES_Optimized)
            blocksOut = DES_ImgProcess_Opt.Engine.encrypt_data(blocksIn, key64, iv64);
            
            % 4. Konwersja powrotna (uint64 -> uint8)
            cipherBytes = typecast(swapbytes(blocksOut), 'uint8');
            
            % 5. Wizualizacja
            len = length(cipherBytes);
            dim = ceil(sqrt(len));
            cipherImgVis = zeros(dim, dim, 'uint8');
            cipherImgVis(1:len) = cipherBytes;
        end
        
        function recoveredImg = decrypt(cipherBytes, keyHex, ivHex, originalSize)
            % 1. Wektoryzacja: uint8 -> uint64
            blocksIn = swapbytes(typecast(cipherBytes, 'uint64'));
            
            key64 = DES_ImgProcess_Opt.hexToUint64(keyHex);
            iv64  = DES_ImgProcess_Opt.hexToUint64(ivHex);
            
            % 2. Szybkie Deszyfrowanie
            blocksOut = DES_ImgProcess_Opt.Engine.decrypt_data(blocksIn, key64, iv64);
            
            % 3. uint64 -> uint8
            decryptedBytes = typecast(swapbytes(blocksOut), 'uint8');
            
            % 4. Unpad
            try
                padLen = double(decryptedBytes(end));
                if padLen > 8 || padLen == 0, error('Bad pad'); end
                decryptedFlat = decryptedBytes(1:end-padLen);
                recoveredImg = reshape(decryptedFlat, originalSize);
            catch
                warning('DES Padding Error');
                recoveredImg = zeros(originalSize, 'uint8');
            end
        end
        
        function val = hexToUint64(hexStr)
            % Konwersja hex string do uint64
            % hex2dec w MATLAB traci precyzję dla > 52 bitów (double).
            % Musimy użyć sscanf lub rozdzielić na części.
            if length(hexStr) > 16, hexStr = hexStr(1:16); end
            
            % Rozbicie na dwie części po 32 bity, aby uniknąć problemów z double
            upper = hex2dec(hexStr(1:8));
            lower = hex2dec(hexStr(9:16));
            
            val = bitor(bitshift(uint64(upper), 32), uint64(lower));
        end
        
        function hexStr = generateRandomHex(byteLength)
             % Generuje losowy klucz hex
             randBytes = randi([0 255], 1, byteLength);
             hexStr = lower(dec2hex(randBytes)');
             hexStr = hexStr(:)';
        end
    end
end