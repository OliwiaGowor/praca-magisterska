classdef ImageBlowfish < handle
    % IMAGEBLOWFISH - Wrapper for Image Encryption
    % Supports: Grayscale, RGB, and Double/Logical inputs.
    
    properties (Access = private)
        cryptoEngine
    end
    
    methods
        function obj = ImageBlowfish(key, pi_hex_str)
            if nargin < 2
                obj.cryptoEngine = Blowfish(key);
            else
                obj.cryptoEngine = Blowfish(key, pi_hex_str);
            end
        end
        
        function [ciphertext, originalSize, originalClass] = encryptImage(obj, imageInput, iv)
            % Input: Filename OR Image Matrix (uint8, double, logical, etc.)
            % Output: 
            %   ciphertext: Encrypted bytes
            %   originalSize: Dimensions [H W] or [H W 3]
            %   originalClass: Class of input (e.g., 'uint8', 'double')
            
            % 1. Load Image
            if ischar(imageInput) || isstring(imageInput)
                img = imread(imageInput);
            else
                img = imageInput;
            end
            
            % 2. Store Metadata
            originalSize = size(img);
            originalClass = class(img);
            
            % 3. Standardize to uint8 (0-255)
            % Matlab 'double' images are 0.0-1.0; 'uint8' are 0-255.
            if isa(img, 'double') || isa(img, 'single')
                img = im2uint8(img); % Handles scaling automatically
            elseif isa(img, 'logical')
                img = uint8(img) * 255;
            elseif ~isa(img, 'uint8')
                % Fallback for uint16, etc.
                img = im2uint8(img);
            end
            
            % 4. Encrypt
            % Flattening is handled inside encryptCBC via (:)'
            ciphertext = obj.cryptoEngine.encryptCBC(img, iv);
        end
        
        function img = decryptImage(obj, ciphertext, iv, originalSize, originalClass)
            % 1. Decrypt
            flatBytes = obj.cryptoEngine.decryptCBC(ciphertext, iv);
            
            % 2. Reshape to Original Dimensions (Works for RGB too)
            try
                img = reshape(flatBytes, originalSize);
            catch
                warning('Decrypted data size (%d) does not match original dimensions.', length(flatBytes));
                img = flatBytes;
                return;
            end
            
            % 3. Restore Original Class
            if nargin > 4
                switch originalClass
                    case 'double'
                        img = im2double(img);
                    case 'logical'
                        img = logical(img > 127);
                    case 'uint16'
                        img = im2uint16(img);
                    % Default is uint8, no action needed
                end
            end
        end
    end
end
