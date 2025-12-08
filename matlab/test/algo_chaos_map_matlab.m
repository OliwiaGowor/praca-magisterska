function [t_enc, t_dec, C1, C2, PT] = algo_chaotic_matlab(data)
    % ALGO_CHAOTIC_MATLAB Benchmark dla 3D Logistic + Improved Chirikov Map
    % Implementacja na podstawie dostarczonego pliku main.m
    
    % 1. Pobranie danych wejściowych
    img_orig = data.img_orig; 
    img_mod  = data.img_mod;
    sz       = data.size;
    [rows, cols, chans] = size(img_orig);
    num_pixels = rows * cols;
    
    % 2. Generowanie kluczy i sekwencji chaotycznej (3D Logistic Map)
    % Parametry z main.m (lub losowe w zakresie stabilności)
    % Klucz dla Chirikova (x, y, k, h)
    key_chirikov = rand(1, 4) * 10; 
    
    % 3. Benchmark SZYFROWANIA
    tic;
    % Faza 1: Generowanie sekwencji i Dyfuzja (Logistic Map)
    [diffused_img, chaotic_seqs] = logistic_diffusion_encrypt(img_orig, num_pixels);
    
    % Faza 2: Konfuzja (Chirikov Map - funkcja zewnętrzna)
    [ct1, ~] = encrypt(diffused_img, key_chirikov);
    t_enc = toc;
    
    % Szyfrowanie obrazu zmodyfikowanego (dla metryk)
    [diffused_mod, ~] = logistic_diffusion_encrypt(img_mod, num_pixels);
    [ct2, ~] = encrypt(diffused_mod, key_chirikov);
    
    % 4. Benchmark DESZYFROWANIA
    tic;
    % Faza 1: Odwrócenie Konfuzji (Inverse Chirikov)
    dec_chirikov = decrypt(ct1, key_chirikov);
    
    % Faza 2: Odwrócenie Dyfuzji (Inverse Logistic Map)
    PT = logistic_diffusion_decrypt(dec_chirikov, chaotic_seqs, sz);
    t_dec = toc;
    
    % Formatowanie wyjścia
    C1 = ct1;
    C2 = ct2;
end

%% --- Funkcje pomocnicze (logika wyciągnięta z main.m) ---

function [out_img, seqs] = logistic_diffusion_encrypt(img, N)
    % Generuje sekwencje i nakłada dyfuzję zgodnie z main.m
    
    % Parametry startowe (z main.m)
    x = zeros(1, N+1); y = zeros(1, N+1); z = zeros(1, N+1);
    x(1)=0.2350; y(1)=0.3500; z(1)=0.7350;
    a=0.0125; b=0.0157; l=3.7700;
    
    % Generowanie ciągu chaotycznego (najwolniejsza część w MATLABie)
    % Optymalizacja: pętla jest konieczna ze względu na zależność rekurencyjną
    for i=1:N
        x(i+1) = l*x(i)*(1-x(i)) + b*y(i)*y(i)*x(i) + a*z(i)*z(i)*z(i);
        y(i+1) = l*y(i)*(1-y(i)) + b*z(i)*z(i)*y(i) + a*x(i)*x(i)*x(i);
        z(i+1) = l*z(i)*(1-z(i)) + b*x(i)*x(i)*z(i) + a*y(i)*y(i)*y(i);
    end
    
    % Usunięcie pierwszego elementu i przycięcie
    x = x(2:end); y = y(2:end); z = z(2:end);
    
    % Zapisz sekwencje do deszyfrowania
    seqs.x = x; seqs.y = y; seqs.z = z;
    
    % Maski całkowitoliczbowe
    Sx = uint8(ceil(mod((x*1000000), 256)));
    Sy = uint8(ceil(mod((y*1000000), 256)));
    Sz = uint8(ceil(mod((z*1000000), 256)));
    
    seqs.Sx = Sx; seqs.Sy = Sy; seqs.Sz = Sz;

    % Operacje na pikselach (z main.m: mnożenie -> uint8 -> XOR)
    % Uwaga: Oryginalny kod wykonuje konwersję uint8(double * uint8), co powoduje utratę danych.
    % Implementujemy dokładnie tak jak w źródle (Lossy process in source logic).
    
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
        % Obsługa grayscale (używamy kanału R/x)
        PR = double(reshape(img(:,:,1), 1, []));
        CDR = x .* PR;
        CCR = bitxor(Sx, uint8(CDR));
        out_img = reshape(CCR, r, c);
    end
end

function out_img = logistic_diffusion_decrypt(img, seqs, sz)
    % Odwraca proces dyfuzji
    [r, c, ch] = size(img);
    num_pixels = r*c;
    
    % Odwracanie sekwencji mnożenia (1/x)
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
