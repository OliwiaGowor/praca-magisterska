classdef DES_ImgProcess_Opt
    % Wrapper handling image conversion using Vectorization
    % and invoking DES_Optimized
    
    properties (Constant)
        Engine = DES_Optimized();
    end
    
    methods (Static)
        function [cipherBytes, cipherImgVis] = encrypt(imgUint8, keyHex, ivHex)
            % 1. Input Vectorization
            % a. PKCS#7 Padding (operations on uint8 vector)
            flatPixels = imgUint8(:)'; 
            padLen = 8 - mod(length(flatPixels), 8);
            paddedPixels = [flatPixels, repmat(uint8(padLen), 1, padLen)];
            
            % b. Vectorization: uint8 -> uint64 (Typecast)
            % MATLAB on PC is usually Little Endian, while DES requires Big Endian.
            % We use swapbytes to correct this.
            blocksIn = swapbytes(typecast(paddedPixels, 'uint64'));
            
            % 2. Prepare keys (hex -> uint64 conversion)
            key64 = DES_ImgProcess_Opt.hexToUint64(keyHex);
            iv64  = DES_ImgProcess_Opt.hexToUint64(ivHex);
            
            % 3. Fast Encryption (DES_Optimized)
            blocksOut = DES_ImgProcess_Opt.Engine.encrypt_data(blocksIn, key64, iv64);
            
            % 4. Convert back (uint64 -> uint8)
            cipherBytes = typecast(swapbytes(blocksOut), 'uint8');
            
            % 5. Visualization
            len = length(cipherBytes);
            dim = ceil(sqrt(len));
            cipherImgVis = zeros(dim, dim, 'uint8');
            cipherImgVis(1:len) = cipherBytes;
        end
        
        function recoveredImg = decrypt(cipherBytes, keyHex, ivHex, originalSize)
            % 1. Vectorization: uint8 -> uint64
            blocksIn = swapbytes(typecast(cipherBytes, 'uint64'));
            
            key64 = DES_ImgProcess_Opt.hexToUint64(keyHex);
            iv64  = DES_ImgProcess_Opt.hexToUint64(ivHex);
            
            % 2. Fast Decryption
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
            % Convert hex string to uint64
            % hex2dec in MATLAB loses precision for > 52 bits (double).
            % We must use sscanf or split into parts.
            if length(hexStr) > 16, hexStr = hexStr(1:16); end
            
            % Split into two 32-bit parts to avoid double precision issues
            upper = hex2dec(hexStr(1:8));
            lower = hex2dec(hexStr(9:16));
            
            val = bitor(bitshift(uint64(upper), 32), uint64(lower));
        end
        
        function hexStr = generateRandomHex(byteLength)
             % Generates a random hex key
             randBytes = randi([0 255], 1, byteLength);
             hexStr = lower(dec2hex(randBytes)');
             hexStr = hexStr(:)';
        end
    end
end
