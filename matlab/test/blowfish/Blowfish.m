classdef Blowfish < handle
    % BLOWFISH - High-Performance Matlab Implementation
    % Fixes:
    %   1. Loads constants from 'BlowfishConstants.mat'
    %   2. Fixes concatenation dimension error in decryptCBC
    
    properties (Access = private)
        P % 18 x 1 uint32
        S % 4 x 256 uint32
        isLittleEndianSystem
    end
    
    methods
        function obj = Blowfish(key)
            % 1. Cache Endianness
            [~, ~, endian] = computer;
            obj.isLittleEndianSystem = isequal(endian, 'L');
            
            % 2. Load Constants
            if exist('BlowfishConstants.mat', 'file')
                data = load('BlowfishConstants.mat');
                obj.P = data.P_init;
                obj.S = data.S_init;
            else
                % Fallback/Error if file missing
                error('BlowfishConstants.mat missing. Please run generate_constants.m first.');
            end
            
            % 3. Key Expansion
            if ischar(key), key = uint8(key); end
            obj.expandKey(key);
        end
        
        function ciphertext = encryptCBC(obj, plaintext, iv)
            if length(iv) ~= 8, error('IV must be 8 bytes'); end
            
            % Force input to row vector for padding
            padded = obj.pkcs7pad(plaintext(:)');
            
            blocks = typecast(padded, 'uint32');
            if obj.isLittleEndianSystem, blocks = swapbytes(blocks); end
            
            % Force blocks to column for consistent processing
            blocks = blocks(:);
            
            ivWords = typecast(uint8(iv), 'uint32');
            if obj.isLittleEndianSystem, ivWords = swapbytes(ivWords); end
            
            L_prev = ivWords(1); 
            R_prev = ivWords(2);
            
            numBlocks = length(blocks) / 2;
            cipherBlocks = zeros(size(blocks), 'uint32');
            
            idx = 1;
            for i = 1:numBlocks
                % blocks is now guaranteed to be column vector or linear
                L_in = bitxor(blocks(idx), L_prev); 
                R_in = bitxor(blocks(idx+1), R_prev);
                
                [L_out, R_out] = obj.encryptBlock(L_in, R_in);
                
                cipherBlocks(idx) = L_out; 
                cipherBlocks(idx+1) = R_out;
                
                L_prev = L_out; 
                R_prev = R_out;
                idx = idx + 2;
            end
            
            if obj.isLittleEndianSystem, cipherBlocks = swapbytes(cipherBlocks); end
            ciphertext = typecast(cipherBlocks, 'uint8');
        end

        function plaintext = decryptCBC(obj, ciphertext, iv)
            blocks = typecast(uint8(ciphertext), 'uint32');
            if obj.isLittleEndianSystem, blocks = swapbytes(blocks); end
            
            % --- FORCE COLUMN VECTOR ---
            blocks = blocks(:); 
            % -----------------------------------------
            
            % Reshape to 2 rows (L top, R bottom)
            blocksMat = reshape(blocks, 2, []);
            
            % Vectorized Block Decryption
            [L_dec, R_dec] = obj.decryptBlock(blocksMat(1,:)', blocksMat(2,:)');
            
            % Prepare IV and Previous Blocks
            ivWords = typecast(uint8(iv), 'uint32');
            if obj.isLittleEndianSystem, ivWords = swapbytes(ivWords); end
            
            % Now this concatenation works because both are columns
            prev = [ivWords(:); blocks(1:end-2)];
            prevMat = reshape(prev, 2, []);
            
            % XOR Decrypted with Previous
            L_plain = bitxor(L_dec, prevMat(1,:)'); 
            R_plain = bitxor(R_dec, prevMat(2,:)');
            
            plainBlocks = [L_plain'; R_plain'];
            
            if obj.isLittleEndianSystem, plainBlocks = swapbytes(plainBlocks(:)); end
            plaintext = obj.pkcs7unpad(typecast(plainBlocks, 'uint8'));
        end
        
        function output = cryptCTR(obj, data, iv)
            if isempty(data), output=[]; return; end
            
            % Force data to column for consistent XOR later
            data = data(:);
            
            iv64 = typecast(uint8(iv), 'uint64');
            if obj.isLittleEndianSystem, iv64 = swapbytes(iv64); end
            
            ctr = iv64 + uint64(0:ceil(length(data)/8)-1)';
            if obj.isLittleEndianSystem, ctr = swapbytes(ctr); end
            
            ctr32 = typecast(ctr, 'uint32');
            ctrMat = reshape(ctr32, 2, []);
            [L, R] = obj.encryptBlock(ctrMat(1,:)', ctrMat(2,:)');
            
            strm = [L'; R'];
            if obj.isLittleEndianSystem, strm = swapbytes(strm(:)); end
            strmBytes = typecast(strm, 'uint8');
            
            % XOR (data is column, strmBytes is column)
            output = bitxor(data, strmBytes(1:length(data)));
            output = output'; % Return as row
        end
    end

    methods (Access = private)
        function [L, R] = encryptBlock(obj, L, R)
            for i = 1:16
                L = bitxor(L, obj.P(i)); R = bitxor(R, obj.F(L));
                tmp = L; L = R; R = tmp;
            end
            tmp = L; L = R; R = tmp;
            R = bitxor(R, obj.P(17)); L = bitxor(L, obj.P(18));
        end

        function [L, R] = decryptBlock(obj, L, R)
            for i = 18:-1:3
                L = bitxor(L, obj.P(i)); R = bitxor(R, obj.F(L));
                tmp = L; L = R; R = tmp;
            end
            tmp = L; L = R; R = tmp;
            R = bitxor(R, obj.P(2)); L = bitxor(L, obj.P(1));
        end

        function y = F(obj, x)
            d = bitand(x, 255); x = bitshift(x, -8);
            c = bitand(x, 255); x = bitshift(x, -8);
            b = bitand(x, 255); x = bitshift(x, -8);
            a = bitand(x, 255);
            s1 = obj.S(1, double(a)+1); s2 = obj.S(2, double(b)+1);
            s3 = obj.S(3, double(c)+1); s4 = obj.S(4, double(d)+1);
            y = bitxor(uint32(mod(double(s1(:))+double(s2(:)), 4294967296)), s3(:));
            y = uint32(mod(double(y)+double(s4(:)), 4294967296));
        end

        function expandKey(obj, key)
            len = length(key);
            if len == 0 || len > 56, error('Key len 1-56'); end
            kIdx = 1;
            for i = 1:18
                val = uint32(0);
                for b = 1:4
                    val = bitor(bitshift(val, 8), uint32(key(mod(kIdx-1, len)+1)));
                    kIdx = kIdx + 1;
                end
                obj.P(i) = bitxor(obj.P(i), val);
            end
            L = uint32(0); R = uint32(0);
            for i = 1:2:18, [L, R] = obj.encryptBlock(L, R); obj.P(i)=L; obj.P(i+1)=R; end
            for i = 1:4, for j = 1:2:256, [L, R] = obj.encryptBlock(L, R); obj.S(i,j)=L; obj.S(i,j+1)=R; end, end
        end
        
        function p = pkcs7pad(~, d), pl=8-mod(length(d),8); if pl==0, pl=8; end; p=[uint8(d(:)'), repmat(uint8(pl),1,pl)]; end
        function d = pkcs7unpad(~, p), pl=double(p(end)); if pl>8||pl==0, error('Bad Pad'); end; d=p(1:end-pl); end
    end
end
