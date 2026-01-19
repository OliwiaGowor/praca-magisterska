function generate_constants()
    % GENERATE_CONSTANTS
    % Calculates standard Blowfish P-array and S-boxes from Pi using BBP.
    % Saves the result to 'BlowfishConstants.mat'.
    
    fprintf('Calculating 8336 hex digits of Pi (this may take 5-10 seconds)...\n');
    
    % 1. Constants needed
    % P: 18 words, S: 4x256=1024 words. Total = 1042 words.
    % Each word is 8 hex digits. Total hex digits = 8336.
    totalHex = 1042 * 8; 
    
    piHex = repmat(' ', 1, totalHex);
    
    % 2. Generate digits
    % Use parallel pool if available for speed, otherwise standard loop
    % (Standard loop shown here for compatibility)
    for k = 1:totalHex
        if mod(k, 1000) == 0, fprintf('...generated %d digits\n', k); end
        val = bbp_pi(k);
        piHex(k) = dec2hex(floor(val * 16));
    end
    
    fprintf('Parsing constants...\n');
    
    % 3. Parse P-array (First 18 words)
    P_init = zeros(18, 1, 'uint32');
    idx = 1;
    for i = 1:18
        chunk = piHex(idx : idx+7);
        P_init(i) = uint32(hex2dec(chunk));
        idx = idx + 8;
    end
    
    % 4. Parse S-boxes (Next 1024 words)
    S_init = zeros(4, 256, 'uint32');
    for i = 1:4
        for j = 1:256
            chunk = piHex(idx : idx+7);
            S_init(i, j) = uint32(hex2dec(chunk));
            idx = idx + 8;
        end
    end
    
    % 5. Save to file
    save('BlowfishConstants.mat', 'P_init', 'S_init');
    fprintf('Success! Saved to BlowfishConstants.mat\n');
end

% --- BBP Helper Functions ---

function digitVal = bbp_pi(n)
    s1 = bbp_sum(1, n-1); s4 = bbp_sum(4, n-1);
    s5 = bbp_sum(5, n-1); s6 = bbp_sum(6, n-1);
    
    frac = 4*s1 - 2*s4 - s5 - s6;
    frac = frac - floor(frac);
    if frac < 0, frac = frac + 1; end
    digitVal = frac;
end

function s = bbp_sum(j, n)
    s = 0;
    % Left Sum
    for k = 0:n
        denom = 8*k + j;
        term = modPow(16, n-k, denom) / denom;
        s = s + term;
        s = s - floor(s);
    end
    % Right Sum
    k = n + 1;
    while true
        denom = 8*k + j;
        term = (16 ^ (n-k)) / denom;
        if term < 1e-16, break; end
        s = s + term;
        s = s - floor(s);
        k = k + 1;
    end
end

function res = modPow(base, exp, modVal)
    res = 1; base = mod(base, modVal);
    while exp > 0
        if mod(exp, 2) == 1, res = mod(res * base, modVal); end
        base = mod(base * base, modVal);
        exp = floor(exp / 2);
    end
end