function [t_enc, t_dec, C1, C2, PT] = algo_aes_matlab(data)
    % AES-256-CBC (Pure MATLAB implementation) - Debug Version
    
    try
        % 1. Generate Key and IV (256-bit key)
        key = uint8(randi([0 255], 1, 32)); 
        iv  = uint8(randi([0 255], 1, 16)); 
        
        % 2. Prepare Data (Force row vector!)
        input_flat     = data.img_flat(:)';      
        input_mod_flat = data.img_flat_mod(:)';
        sz             = data.size;
        
        % 3. ENCRYPTION
        tic;
        ct1_flat = aes_mode(input_flat, key, iv, 'encrypt', 'cbc');
        t_enc = toc;
        
        ct2_flat = aes_mode(input_mod_flat, key, iv, 'encrypt', 'cbc');
        
        % 4. DECRYPTION
        tic;
        pt_flat_padded = aes_mode(ct1_flat, key, iv, 'decrypt', 'cbc');
        t_dec = toc;
        
        % 5. Format Results
        num_pixels = prod(sz);
        % Crop to original size (remove padding if present) for reshaping
        C1 = reshape(ct1_flat(1:min(end,num_pixels)), sz);
        C2 = reshape(ct2_flat(1:min(end,num_pixels)), sz);
        PT = reshape(pt_flat_padded(1:min(end,num_pixels)), sz);
        
    catch ME
        % Print full error stack to see exact problem location
        fprintf('Error in algo_aes_matlab details:\n');
        fprintf('%s\n', getReport(ME));
        rethrow(ME); % Pass error up to run_benchmarks
    end
end
