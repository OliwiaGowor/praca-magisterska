function [t_enc, t_dec, C1, C2, PT] = algo_chaos_map_matlab(data)
    % ALGO_CHAOS_MAP_MATLAB Benchmark for 3D Logistic + Improved Chirikov Map
    % Implementation based on the provided main logic.
    
    % 1. Input data retrieval
    img_orig = data.img_orig; 
    img_mod  = data.img_mod;
    sz       = data.size;
    [rows, cols, chans] = size(img_orig);
    num_pixels = rows * cols;
    
    % 2. Key generation and chaotic sequence (3D Logistic Map)
    % Parameters from main.m (or random within stability range)
    % Key for Chirikov (x, y, k, h)
    key_chirikov = rand(1, 4) * 10; 
    
    % 3. ENCRYPTION Benchmark
    tic;
    % Phase 1: Sequence Generation and Diffusion (Logistic Map)
    [diffused_img, chaotic_seqs] = logistic_diffusion_encrypt(img_orig, num_pixels);
    
    % Phase 2: Confusion (Chirikov Map - external function)
    [ct1, ~] = encrypt(diffused_img, key_chirikov);
    t_enc = toc;
    
    % Encrypt modified image (for metrics)
    [diffused_mod, ~] = logistic_diffusion_encrypt(img_mod, num_pixels);
    [ct2, ~] = encrypt(diffused_mod, key_chirikov);
    
    % 4. DECRYPTION Benchmark
    tic;
    % Phase 1: Inverse Confusion (Inverse Chirikov)
    dec_chirikov = decrypt(ct1, key_chirikov);
    
    % Phase 2: Inverse Diffusion (Inverse Logistic Map)
    PT = logistic_diffusion_decrypt(dec_chirikov, chaotic_seqs, sz);
    t_dec = toc;
    
    % Output formatting
    C1 = ct1;
    C2 = ct2;
end

%% --- Helper Functions (Logic extracted from main.m) ---

function [out_img, seqs] = logistic_diffusion_encrypt(img, N)
    % Generates sequences and applies diffusion according to main.m
    
    % Start parameters (from main.m)
    x = zeros(1, N+1); y = zeros(1, N+1); z = zeros(1, N+1);
    x(1)=0.2350; y(1)=0.3500; z(1)=0.7350;
    a=0.0125; b=0.0157; l=3.7700;
    
    % Generate chaotic sequence (Slowest part in MATLAB)
    % Optimization: Loop is necessary due to recursive dependency
    for i=1:N
        x(i+1) = l*x(i)*(1-x(i)) + b*y(i)*y(i)*x(i) + a*z(i)*z(i)*z(i);
        y(i+1) = l*y(i)*(1-y(i)) + b*z(i)*z(i)*y(i) + a*x(i)*x(i)*x(i);
        z(i+1) = l*z(i)*(1-z(i)) + b*x(i)*x(i)*z(i) + a*y(i)*y(i)*y(i);
    end
    
    % Remove first element and trim
    x = x(2:end); y = y(2:end); z = z(2:end);
    
    % Save sequences for decryption
    seqs.x = x; seqs.y = y; seqs.z = z;
    
    % Integer masks
    Sx = uint8(ceil(mod((x*1000000), 256)));
    Sy = uint8(ceil(mod((y*1000000), 256)));
    Sz = uint8(ceil(mod((z*1000000), 256)));
    
    seqs.Sx = Sx; seqs.Sy = Sy; seqs.Sz = Sz;

    % Pixel operations (from main.m: multiply -> uint8 -> XOR)
    % Note: Original code performs uint8(double * uint8) conversion which causes data loss.
    % We implement exactly as in the source logic.
    
    [r, c, ch] = size(img);
    if ch >= 3
        PR = double(reshape(img(:,:,1), 1, []));
        PG = double(reshape(img(:,:,2), 1, []));
        PB = double(reshape(img(:,:,3), 1, []));
        
        CDR = x .* PR; 
        CDG = y .* PG;
        CDB = z .* PB;
        
        CCR = bitxor(Sx, uint8(CDR));
        CCG = bitxor(Sy, uint8(CDG));
        CCB = bitxor(Sz, uint8(CDB));
        
        out_img = cat(3, reshape(CCR,r,c), reshape(CCG,r,c), reshape(CCB,r,c));
    else
        % Grayscale handling (using R/x channel)
        PR = double(reshape(img(:,:,1), 1, []));
        CDR = x .* PR;
        CCR = bitxor(Sx, uint8(CDR));
        out_img = reshape(CCR, r, c);
    end
end

function out_img = logistic_diffusion_decrypt(img, seqs, sz)
    % Reverses the diffusion process
    [r, c, ch] = size(img);
    num_pixels = r*c;
    
    % Reverse multiplication sequence (1/x)
    xi = 1./seqs.x; yi = 1./seqs.y; zi = 1./seqs.z;
    
    if ch >= 3
        DR = reshape(img(:,:,1), 1, []);
        DG = reshape(img(:,:,2), 1, []);
        DB = reshape(img(:,:,3), 1, []);
        
        % XOR (Odwrotność XOR to XOR)
        DDR = bitxor(seqs.Sx, uint8(DR));
        DDG = bitxor(seqs.Sy, uint8(DG));
        DDB = bitxor(seqs.Sz, uint8(DB));
        
        % Odwrócenie mnożenia
        DDDR = xi .* double(DDR);
        DDDG = yi .* double(DDG);
        DDDB = zi .* double(DDB);
        
        out_img = uint8(cat(3, reshape(DDDR,r,c), reshape(DDDG,r,c), reshape(DDDB,r,c)));
    else
        DR = reshape(img(:,:,1), 1, []);
        DDR = bitxor(seqs.Sx, uint8(DR));
        DDDR = xi .* double(DDR);
        out_img = uint8(reshape(DDDR, r, c));
    end
end
