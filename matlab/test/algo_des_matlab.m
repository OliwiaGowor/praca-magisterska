function [t_enc, t_dec, C1, C2, PT] = algo_des_matlab(data)
    % DES-CBC (Pure MATLAB)
    
    % Generate keys
    key_hex = DES_ImgProcess_Opt.generateRandomHex(8);
    iv_hex  = DES_ImgProcess_Opt.generateRandomHex(8);
    
    % Extract data from structure to simple local variables.
    local_img_orig = data.img_orig; 
    local_img_mod  = data.img_mod;
    local_size     = data.size;
    
    tic;
    % Use the local variable for encryption
    [ct_bytes_1, C1] = DES_ImgProcess_Opt.encrypt(local_img_orig, key_hex, iv_hex);
    t_enc = toc;
    
    % Encrypt modified image
    [~, C2] = DES_ImgProcess_Opt.encrypt(local_img_mod, key_hex, iv_hex);
    
    tic;
    % Decryption
    PT = DES_ImgProcess_Opt.decrypt(ct_bytes_1, key_hex, iv_hex, local_size);
    t_dec = toc;
end
