function C_attacked = apply_attack(C_orig, attack_type, intensity, img_tamper)
    % APPLY_ATTACK Simulates alteration attacks on ciphertext.
    % Supports: crop, contiguous_crop, tamper, salt_pepper, gaussian, speckle, poisson, bit_crop.
    
        C_attacked = C_orig;
        [rows, cols] = size(C_orig);
        num_pixels = rows * cols;
        
        % Convert to uint8/double as needed (noise functions require compatibility)
        if isa(C_attacked, 'uint8')
            C_attacked_u8 = C_attacked;
        else
            C_attacked_u8 = im2uint8(C_attacked);
        end
    
        switch attack_type
            case 'crop'
                % Scattered attack (Random Cropping) -> Noise on ciphertext
                num_attacked_pixels = round(num_pixels * intensity);
                idx_to_crop = randperm(num_pixels, num_attacked_pixels);
                C_attacked(idx_to_crop) = 0; 
    
            case 'contiguous_crop'
                % Contiguous attack (Block Cropping) -> Black hole on ciphertext
                
                % Calculate side length of the square to crop
                side = round(sqrt(intensity * num_pixels));
                
                % Image center
                r_start = round(rows/2 - side/2);
                c_start = round(cols/2 - side/2);
                
                % Boundary check
                if r_start < 1, r_start = 1; end
                if c_start < 1, c_start = 1; end
                r_end = r_start + side;
                c_end = c_start + side;
                if r_end > rows, r_end = rows; end
                if c_end > cols, c_end = cols; end
                
                % If data is a vector (1D)
                if isvector(C_attacked)
                    center_idx = round(numel(C_attacked)/2);
                    len = round(numel(C_attacked) * intensity);
                    start_idx = max(1, center_idx - floor(len/2));
                    end_idx = min(numel(C_attacked), center_idx + ceil(len/2));
                    C_attacked(start_idx:end_idx) = 0;
                else
                    % For 2D matrix
                    C_attacked(r_start:r_end, c_start:c_end) = 0;
                end
    
            case 'bit_crop'
                % Zeroing least significant bits (or random planes)
                bits_to_remove = round(intensity * 8);
                if bits_to_remove < 1, bits_to_remove = 1; end
                planes = randperm(8, bits_to_remove);
                for p = planes
                    C_attacked = bitset(C_attacked, p, 0);
                end
    
            case 'tamper'
                % Substituting a fragment with another image (Tampering)
                if nargin < 4 || isempty(img_tamper)
                    img_tamper = uint8(randi([0 255], rows, cols));
                end
                if size(img_tamper,1) ~= rows || size(img_tamper,2) ~= cols
                    img_tamper = imresize(img_tamper, [rows, cols]);
                end
                if size(img_tamper,3) > 1
                    img_tamper = rgb2gray(img_tamper);
                end
    
                num_attacked_pixels = round(num_pixels * intensity);
                idx_to_tamper = randperm(num_pixels, num_attacked_pixels);
                C_attacked(idx_to_tamper) = img_tamper(idx_to_tamper);
    
            case 'salt_pepper'
                C_attacked = imnoise(C_attacked_u8, 'salt & pepper', intensity);
            case 'gaussian'
                C_attacked = imnoise(C_attacked_u8, 'gaussian', 0, intensity);
            case 'speckle'
                C_attacked = imnoise(C_attacked_u8, 'speckle', intensity);
            case 'poisson'
                C_attacked = imnoise(C_attacked_u8, 'poisson');
                
            otherwise
                error('Nieznany typ ataku: %s', attack_type);
        end
        
        % Restore type
        if isa(C_orig, 'uint8')
            C_attacked = uint8(C_attacked);
        else
            C_attacked = double(C_attacked);
        end
    end