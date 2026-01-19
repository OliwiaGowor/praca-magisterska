function [out1, out2] = chaos_2d_adapter(mode, varargin)
    % CHAOS_2D_ADAPTER
    % Exact cryptographic core of Algorithms 3 & 4
    % Reference: Entropy 2024, 26, 603
    
        switch lower(mode)
            case 'encrypt'
                [out1, out2] = encrypt_core(varargin{1});
            case 'decrypt'
                out1 = decrypt_core(varargin{1}, varargin{2});
                out2 = [];
            otherwise
                error('Mode must be "encrypt" or "decrypt"');
        end
    end
    
    % =========================================================
    % ENCRYPT (Algorithm 3)
    % =========================================================
    function [C, keys] = encrypt_core(P)
    
    params.a  = 20;
    params.b  = 30;
    params.x0 = 0.123987;
    params.y0 = 0.987321;
    params.Pre = uint8(50);
    
    if size(P,3)==3, P = rgb2gray(P); end
    P = uint8(P);
    [M,N] = size(P);
    L = M*N;
    
    x = params.x0; y = params.y0;
    X = zeros(1,L); Y = zeros(1,L);
    for i=1:L
        [x,y] = HyperchaoticMap(x,y,params.a,params.b);
        X(i)=x; Y(i)=y;
    end
    
    [~, perm] = sort(X);       % perm is 1:L
    invp = zeros(1,L);
    invp(perm) = 1:L;
    
    K = uint8(mod(floor(Y*1e14),256));
    
    Pp = P(:);
    Pp = Pp(perm);
    
    C1 = zeros(L,1,'uint8');
    Pre = params.Pre;
    for i=1:L
        C1(i) = bitxor(uint8(mod(double(Pp(i))+double(Pre),256)),K(i));
        Pre = C1(i);
    end
    
    C2 = zeros(L,1,'uint8');
    Pre = params.Pre;
    for i=L:-1:1
        C2(i) = bitxor(uint8(mod(double(C1(i))+double(Pre),256)),K(i));
        Pre = C2(i);
    end
    
    C = reshape(C2,M,N);
    
    keys.K = K;
    keys.perm = perm;
    keys.inv_perm = invp;
    keys.params = params;
    end
    
    % =========================================================
    % DECRYPT (Algorithm 4)
    % =========================================================
    function P = decrypt_core(C, keys)
    
    params = keys.params;
    K      = keys.K;
    perm   = keys.perm;
    invp   = keys.inv_perm;
    
    [M,N] = size(C);
    L = M*N;
    Cvec = C(:);
    
    C1 = zeros(L,1,'uint8');
    Pre = params.Pre;
    for i=L:-1:1
        C1(i) = uint8(mod(double(bitxor(Cvec(i),K(i))) - double(Pre),256));
        Pre = Cvec(i);
    end
    
    Pp = zeros(L,1,'uint8');
    Pre = params.Pre;
    for i=1:L
        Pp(i) = uint8(mod(double(bitxor(C1(i),K(i))) - double(Pre),256));
        Pre = C1(i);
    end
    
    Pvec = Pp(invp);
    P = reshape(Pvec,M,N);
    end
    
    function [x,y] = HyperchaoticMap(x,y,a,b)
        x = sin(a*pi*x + b*y)^2;
        if y==0
            y = 0;
        else
            y = cos(b*pi/y + a*x)^2;
        end
    end
